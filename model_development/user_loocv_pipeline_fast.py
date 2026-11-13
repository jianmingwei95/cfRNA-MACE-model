# ============================================================
# Original "cfrna loocv pipeline.ipynb", run unchanged EXCEPT:
#   - step 1 of select_top_n_features (univariable lifelines Cox C-index for every gene in every fold)
#     is computed once for all folds by a vectorised leave-one-out update; verified against lifelines below.
#   - N_FEATURES taken from the command line; outputs written to 03_results/loocv_user_code/
#   - an additional, separately reported "centred" risk score (held-out score minus the training-fold mean score)
# Model, correlation filter, fallback rule, CoxPHSurvivalAnalysis(alpha=0.1), uncentred predict()
# and concordance_index_censored are identical to the notebook.
# Usage: /opt/miniconda3/bin/python user_loocv_pipeline_fast.py <N_FEATURES>
# ============================================================
import sys, os, numpy as np, pandas as pd, warnings
from collections import Counter
from sklearn.model_selection import LeaveOneOut
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.util import Surv
from sksurv.metrics import concordance_index_censored
from sklearn.metrics import roc_auc_score, average_precision_score
from lifelines import CoxPHFitter
warnings.filterwarnings("ignore")

N_FEATURES = int(sys.argv[1]) if len(sys.argv) > 1 else 5
OUT = "path/to/project/external/03_results/loocv_user_code"
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv("path/to/project/data/survival_data_cleanqc.csv")
X = df.drop(columns=["EP_major", "EP_major_time", "STUDY_NUMBER"])
y = df["EP_major"].astype(int)
time = df["EP_major_time"]
n, p = X.shape
print(f"Samples {n} | events {y.sum()} | features {p} | N_FEATURES {N_FEATURES}")

# ------------------------------------------------------------
# Vectorised leave-one-out univariable C-index (replaces per-fold lifelines loop)
# lifelines' concordance_index_ of a univariable Cox model = Harrell's C of beta*x, i.e. C(x) if beta > 0 else 1 - C(x).
# Pairs: i event, t_j > t_i, or t_j == t_i with j censored. Removing patient k removes only pairs involving k.
# ------------------------------------------------------------
Xv = X.values.astype(float); t = time.values.astype(float); e = y.values.astype(int)
pi, pj = [], []
for i in np.where(e == 1)[0]:
    js = np.where((t > t[i]) | ((t == t[i]) & (e == 0)))[0]
    pi.extend([i] * len(js)); pj.extend(js)
pi = np.array(pi); pj = np.array(pj); D = len(pi)
dk = np.bincount(np.concatenate([pi, pj]), minlength=n)
Ntot = np.zeros(p); Nk = np.zeros((n, p))
for s in range(0, p, 400):
    A = Xv[pi, s:s+400]; B = Xv[pj, s:s+400]
    V = (A > B) + 0.5 * (A == B)
    Ntot[s:s+400] = V.sum(0)
    np.add.at(Nk[:, s:s+400], pi, V); np.add.at(Nk[:, s:s+400], pj, V)
CLOO = (Ntot[None, :] - Nk) / (D - dk)[:, None]          # C(x) without patient k
BETA_POS = None                                         # sign of beta decided per fold below via C side

def loo_cindex_dict(k, train_idx):
    c = CLOO[k]
    # direction from univariable Cox: beta > 0  <=>  C(x) > 0.5 (verified against lifelines below)
    return dict(zip(X.columns, np.where(c >= 0.5, c, 1 - c)))

# ---- verification against the original lifelines computation on one fold ----
rng = np.random.default_rng(7); kk = int(rng.integers(n)); tr0 = np.setdiff1d(np.arange(n), [kk])
genes = list(rng.choice(X.columns, 40, replace=False))
fast = loo_cindex_dict(kk, tr0); diffs = []
for g in genes:
    d0 = pd.DataFrame({"time": t[tr0], "event": e[tr0], "feature": Xv[tr0, X.columns.get_loc(g)]})
    try:
        cph = CoxPHFitter(); cph.fit(d0, duration_col="time", event_col="event"); diffs.append(abs(cph.concordance_index_ - fast[g]))
    except Exception: pass
print(f"Verification fold {kk}: {len(diffs)} genes, max |fast - lifelines| = {max(diffs):.2e}")

def select_top_n_features(X_train, univariate_cindices, n=10, corr_threshold=0.7):
    sorted_features = sorted(univariate_cindices.items(), key=lambda x: x[1], reverse=True)   # identical to notebook
    selected = []
    for feat, c_idx in sorted_features:
        if len(selected) >= n: break
        if len(selected) == 0: selected.append(feat)
        else:
            max_corr = max([abs(X_train[feat].corr(X_train[sf])) for sf in selected])
            if max_corr < corr_threshold: selected.append(feat)
    if len(selected) < n:
        for feat, c_idx in sorted_features:
            if feat not in selected:
                selected.append(feat)
                if len(selected) >= n: break
    return selected[:n]

loo = LeaveOneOut()
risk_scores, risk_centred, y_true, time_true, selected_all = [], [], [], [], []
feature_counter = Counter()
for fold_id, (train_idx, test_idx) in enumerate(loo.split(X), start=1):
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
    time_train, time_test = time.iloc[train_idx], time.iloc[test_idx]
    sel = select_top_n_features(X_train, loo_cindex_dict(test_idx[0], train_idx), n=N_FEATURES, corr_threshold=0.7)
    selected_all.append(sel); feature_counter.update(sel)
    y_surv_train = Surv.from_arrays(event=y_train.astype(bool).values, time=time_train.values)
    try:
        model = CoxPHSurvivalAnalysis(alpha=0.1)
        model.fit(X_train[sel].values, y_surv_train)
        rs = model.predict(X_test[sel].values)[0]
        rc = rs - model.predict(X_train[sel].values).mean()
    except Exception:
        rs = rc = np.nan
    risk_scores.append(rs); risk_centred.append(rc)
    y_true.append(y_test.values[0]); time_true.append(time_test.values[0])
    if fold_id % 50 == 0: print(f"  fold {fold_id}/{n}", flush=True)

y_true = np.array(y_true); time_true = np.array(time_true)
res = {}
for lab, rsk in [("original (uncentred, as notebook)", np.array(risk_scores)), ("centred (training-fold mean removed)", np.array(risk_centred))]:
    v = ~np.isnan(rsk)
    c = concordance_index_censored(y_true[v].astype(bool), time_true[v], rsk[v])[0]
    rng2 = np.random.default_rng(1); bs = []
    idx = np.where(v)[0]
    for _ in range(1000):
        b = rng2.choice(idx, len(idx), replace=True)
        if y_true[b].sum() > 0: bs.append(concordance_index_censored(y_true[b].astype(bool), time_true[b], rsk[b])[0])
    res[lab] = (c, np.percentile(bs, 2.5), np.percentile(bs, 97.5), roc_auc_score(y_true[v], rsk[v]), v.sum())
    print(f"[{lab}] N_FEATURES={N_FEATURES}: C-index {c:.3f} (95% CI {np.percentile(bs,2.5):.3f}-{np.percentile(bs,97.5):.3f}) | ROC-AUC {res[lab][3]:.3f} | valid {v.sum()}/{n}")

fs = pd.DataFrame([{"Feature": f, "Count": c, "Frequency": c / n} for f, c in feature_counter.most_common()])
print("Top 12 selected features:"); print(fs.head(12).to_string(index=False))
pd.DataFrame([{"N_Features": N_FEATURES, "Score": k, "C_Index": v[0], "CI_low": v[1], "CI_high": v[2], "ROC_AUC": v[3], "Valid": v[4]} for k, v in res.items()]).to_csv(f"{OUT}/user_loocv_{N_FEATURES}features_performance.csv", index=False)
fs.to_csv(f"{OUT}/user_loocv_{N_FEATURES}features_feature_stability.csv", index=False)
pd.DataFrame({"Fold": range(1, n + 1), "Event": y_true, "Time": time_true, "RiskScore": risk_scores, "RiskScore_centred": risk_centred,
              "Features": [",".join(f) for f in selected_all]}).to_csv(f"{OUT}/user_loocv_{N_FEATURES}features_fold_details.csv", index=False)
