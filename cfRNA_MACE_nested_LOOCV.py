# ============================================================
# STRICT NESTED LOOCV
# + Per-Fold Independent Top 10 Feature Selection
# + Completely Leakage-Free
# + Feature Ranking Using C-index
# ============================================================

import numpy as np
import pandas as pd
from collections import Counter
from sklearn.model_selection import LeaveOneOut
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.util import Surv
from sksurv.metrics import concordance_index_censored
from sklearn.metrics import roc_auc_score, average_precision_score
from lifelines import CoxPHFitter
from tqdm import tqdm
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

print("="*80)
print(" 🚀 STRICT NESTED LOOCV: Per-Fold Feature Selection (Top 10)")
print("="*80)

# ============================================================
# Load Data
# ============================================================
df = pd.read_csv("survival_data_cleanqc.csv")
X = df.drop(columns=["EP_major", "EP_major_time", "STUDY_NUMBER"])
y = df["EP_major"].astype(int)
time = df["EP_major_time"]

print(f"\n📊 Dataset:")
print(f"   Samples : {len(X)}")
print(f"   Events  : {y.sum()} ({y.mean()*100:.1f}%)")
print(f"   Features: {X.shape[1]}")

# ============================================================
# Core Function: Select Top N Features on Training Set
# ============================================================
def select_top_n_features(X_train, time_train, y_train, n=10, corr_threshold=0.7):
    """
    Select Top N low-correlated features on the training set
    
    Steps:
    1. Compute univariate C-index for each feature using Cox models
    2. Rank features by descending C-index
    3. Greedy selection prioritising high C-index and low correlation
    4. Relax correlation constraints if fewer than n features remain
    
    Parameters:
        X_train: Training feature matrix
        time_train: Training survival time
        y_train: Training event labels
        n: Number of features to select
        corr_threshold: Correlation threshold
    
    Returns:
        selected: List of selected features
    """
    
    # Step 1: Compute univariate C-index
    univariate_cindices = {}
    
    for col in X_train.columns:
        try:
            df_cox = pd.DataFrame({
                'time': time_train.values,
                'event': y_train.values,
                'feature': X_train[col].values
            })
            
            cph = CoxPHFitter()
            cph.fit(df_cox, duration_col='time', event_col='event')
            
            univariate_cindices[col] = cph.concordance_index_
        except:
            univariate_cindices[col] = 0.5
    
    # Step 2: Rank by C-index
    sorted_features = sorted(
        univariate_cindices.items(),
        key=lambda x: x[1],
        reverse=True
    )
    
    # Step 3: Greedy selection of low-correlated features
    selected = []
    
    for feat, c_idx in sorted_features:
        if len(selected) >= n:
            break
        
        if len(selected) == 0:
            selected.append(feat)
        else:
            max_corr = max([
                abs(X_train[feat].corr(X_train[sf]))
                for sf in selected
            ])
            
            if max_corr < corr_threshold:
                selected.append(feat)
    
    # Step 4: Relax correlation constraint if needed
    if len(selected) < n:
        for feat, c_idx in sorted_features:
            if feat not in selected:
                selected.append(feat)
                if len(selected) >= n:
                    break
    
    return selected[:n]

# ============================================================
# Main Loop: Strict Nested LOOCV
# ============================================================
N_FEATURES = 10

print(f"\n⚙️  Strategy:")
print(f"   ├─ Outer CV: LOOCV ({len(X)} folds)")
print(f"   ├─ Inner: Per-fold feature selection")
print(f"   ├─ Features per fold: {N_FEATURES}")
print(f"   ├─ Selection method: Univariate C-index (Cox)")
print(f"   ├─ Correlation threshold: 0.7")
print(f"   └─ Model: Cox PH with L2 regularisation (alpha=0.1)")

print(f"\n⏳ Running Strict Nested LOOCV...")
print(f"⚠️  Note: This may take ~2–3 hours due to repeated Cox model fitting.")

loo = LeaveOneOut()

risk_scores = []
y_true = []
time_true = []
selected_features_all = []
feature_counter = Counter()

for fold_id, (train_idx, test_idx) in enumerate(
    tqdm(loo.split(X), total=len(X), desc="LOOCV"),
    start=1
):
    X_train = X.iloc[train_idx]
    X_test = X.iloc[test_idx]
    y_train = y.iloc[train_idx]
    y_test = y.iloc[test_idx]
    time_train = time.iloc[train_idx]
    time_test = time.iloc[test_idx]
    
    # ============================================================
    # Step 1: Feature Selection on Training Set Only
    # ============================================================
    selected_features = select_top_n_features(
        X_train, time_train, y_train,
        n=N_FEATURES,
        corr_threshold=0.7
    )
    
    selected_features_all.append(selected_features)
    feature_counter.update(selected_features)
    
    # ============================================================
    # Step 2: Train Model Using Selected Features
    # ============================================================
    X_train_selected = X_train[selected_features]
    X_test_selected = X_test[selected_features]
    
    y_surv_train = Surv.from_arrays(
        event=y_train.astype(bool).values,
        time=time_train.values
    )
    
    try:
        model = CoxPHSurvivalAnalysis(alpha=0.1)
        model.fit(X_train_selected.values, y_surv_train)
        
        risk_score = model.predict(X_test_selected.values)[0]
    except:
        risk_score = np.nan
    
    risk_scores.append(risk_score)
    y_true.append(y_test.values[0])
    time_true.append(time_test.values[0])

y_true = np.array(y_true)
time_true = np.array(time_true)
risk_scores = np.array(risk_scores)

print("\n✅ Nested LOOCV completed!\n")
