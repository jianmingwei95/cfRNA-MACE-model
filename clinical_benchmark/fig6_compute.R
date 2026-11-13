## ============================================================================
## Figure 6 data: clinical utility of the locked cfRNA score, replication cohort (serum, n = 277)
##   All parameters fixed in discovery (n = 344) and applied unchanged to replication:
##   - high-risk group = top tertile of the locked-model score in discovery (pre-specified threshold,
##     locked_threshold.csv, as used for sensitivity/specificity in Figure 4)
##   - SMART2 + cfRNA combination weights (Cox, MACE, discovery; locked-model score); comparator = SMART2
##   - absolute 3-year ischaemic-stroke risk: slope + baseline hazard from discovery Cox models
##   A/B cumulative incidence by risk group (ischaemic stroke / MACE); C time-dependent AUC (IS);
##   D calibration at 3 years (IS); E decision curve at 3 years (IS); F HR per SD by sex, both
##   cohorts and end points, with interaction P.
## Outputs: 03_results/fig6/*.csv (aggregate curves and statistics only)
## ============================================================================
root <- "path/to/project/external"
src <- readLines(file.path(root, "clinical_benchmark", "benchmark_clinical_scores.R"))
eval(parse(text = src[seq_len(grep("^## ---------------- 0\\) cohorts", src) - 1)]))
suppressMessages({ library(glmnet); library(timeROC); library(dcurves) })
OUT6 <- file.path(root, "03_results", "fig6"); sup <- function(x) ifelse(x < 10, "<10", as.character(x))

## ---- inputs --------------------------------------------------------------------------------------
Lk <- readRDS(file.path(root, "03_results", "disc344_to_external", "model_5gene.rds"))
Dg <- fread("path/to/project/data/survival_data_cleanqc.csv", select = c("STUDY_NUMBER", Lk$sig))
Dg[, lock := as.numeric(predict(Lk$fit, as.matrix(Dg[, Lk$sig, with = FALSE]), type = "link"))]
lo <- fread("path/to/project/data/strict_nested_loocv_merged.csv", sep = ";", dec = ",")
old <- fread("path/to/project/data/clinical_imp_final.csv")
disc <- merge(merge(old[, .(STUDY_NUMBER)], Dg[, .(STUDY_NUMBER, lock)], by = "STUDY_NUMBER"), lo[, .(STUDY_NUMBER, loocv = RiskScore)], by = "STUDY_NUMBER")
rep <- fread(file.path(root, "03_results", "disc344_to_external", "external3batch_riskscore_and_genes.csv"))[, .(STUDY_NUMBER, lock = risk)]
disc$SMART <- scores(impute(disc$STUDY_NUMBER)$data)$SMART2      # column name kept generic: holds SMART2
rep$SMART  <- scores(impute(rep$STUDY_NUMBER)$data)$SMART2
sx <- function(ids) as.numeric(zap_labels(cl$Gender[match(ids, cl$STUDY_NUMBER)]))
disc$male <- sx(disc$STUDY_NUMBER); rep$male <- sx(rep$STUDY_NUMBER)
for (ep in c("EP_major", "EP_ischemic_stroke")) { a <- outcome(disc$STUDY_NUMBER, ep); b <- outcome(rep$STUDY_NUMBER, ep)
  disc[[paste0(ep, "_t")]] <- a$t; disc[[paste0(ep, "_e")]] <- a$e; rep[[paste0(ep, "_t")]] <- b$t; rep[[paste0(ep, "_e")]] <- b$e }

## ---- locked parameters from discovery -------------------------------------------------------------
cut <- unname(quantile(disc$lock, c(1/3, 2/3)))
fm <- coxph(Surv(EP_major_t, EP_major_e) ~ SMART + lock, data = disc); wA <- coef(fm)[["SMART"]]; wB <- coef(fm)[["lock"]]   # locked-model cfRNA score in both cohorts
disc <- copy(disc); rep <- copy(rep); disc[, comb := wA * SMART + wB * lock]; rep[, comb := wA * SMART + wB * lock]
fwrite(data.table(cut_low = cut[1], cut_high = cut[2], w_SMART = wA, w_cfRNA = wB), file.path(OUT6, "fig6_locked_parameters.csv"))
cat(sprintf("[locked] tertile cuts %.3f / %.3f | combined = %.3f x SMART2 + %.3f x cfRNA\n", cut[1], cut[2], wA, wB))
thr <- fread(file.path(root, "03_results", "fig4", "locked_threshold.csv"))$threshold_top_tertile
stopifnot(abs(thr - cut[2]) < 1e-6)
grp <- function(x) factor(ifelse(x > thr, "High", "Low"), c("Low", "High"))
rep[, group := grp(lock)]; cat("[groups] replication sizes:", paste(levels(rep$group), sup(as.integer(table(rep$group))), sep = "=", collapse = " "), "\n")

## ---- A/B cumulative incidence by risk group -----------------------------------------------------
KMo <- list(); KS <- list(); RK <- list()
for (ep in c("EP_ischemic_stroke", "EP_major")) {
  tt <- rep[[paste0(ep, "_t")]]; ee <- rep[[paste0(ep, "_e")]]
  f <- survfit(Surv(tt, ee) ~ group, data = rep)
  s <- summary(f, censored = TRUE)
  KMo[[ep]] <- data.table(endpoint = ep, group = sub("group=", "", s$strata), time = s$time, cuminc = 1 - s$surv, lo = 1 - s$upper, hi = 1 - s$lower)
  ar <- summary(f, times = 0:5, extend = TRUE)
  RK[[ep]] <- data.table(endpoint = ep, group = sub("group=", "", ar$strata), time = ar$time, n_risk = ar$n.risk)
  cx <- coxph(Surv(tt, ee) ~ group, data = rep); lr <- survdiff(Surv(tt, ee) ~ group, data = rep)
  ci <- exp(confint(cx)); s5 <- summary(f, times = 5, extend = TRUE)$surv
  KS[[ep]] <- data.table(endpoint = ep, HR_high = exp(coef(cx))[[1]], HR_high_lo = ci[1, 1], HR_high_hi = ci[1, 2],
                         p_cox = summary(cx)$coefficients[, 5], p_logrank = 1 - pchisq(lr$chisq, 1),
                         cuminc5_low = 1 - s5[1], cuminc5_high = 1 - s5[2])
}
fwrite(rbindlist(KMo), file.path(OUT6, "fig6ab_km.csv")); fwrite(rbindlist(RK), file.path(OUT6, "fig6ab_atrisk.csv")); fwrite(rbindlist(KS), file.path(OUT6, "fig6ab_stats.csv"))
cat("[A/B]\n"); print(rbindlist(KS)[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])

## ---- C time-dependent AUC, ischaemic stroke -------------------------------------------------------
TS <- c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
tA <- rbindlist(lapply(c(SMART2 = "SMART", cfRNA = "lock", `SMART2 + cfRNA` = "comb"), function(v) {
  r <- timeROC(T = rep$EP_ischemic_stroke_t, delta = rep$EP_ischemic_stroke_e, marker = rep[[v]], cause = 1, weighting = "marginal", times = TS, iid = TRUE)
  ci <- confint(r, level = 0.95)$CI_AUC / 100
  data.table(time = TS, AUC = r$AUC, lo = ci[, 1], hi = ci[, 2]) }), idcol = "model")
fwrite(tA, file.path(OUT6, "fig6c_timeauc.csv")); cat("[C] time-AUC at 3 y:\n"); print(tA[time == 3, .(model, AUC = round(AUC, 3), lo = round(lo, 3), hi = round(hi, 3))])

## ---- D calibration at 3 years (IS) + E decision curve -------------------------------------------------
H <- 3
## slope fixed from the MACE Cox model in discovery (same model that defines the score); only the 3-year
## ischaemic-stroke baseline hazard is taken from discovery (Breslow, slope held fixed as an offset)
risk3 <- function(v) { f <- coxph(as.formula(paste("Surv(EP_major_t, EP_major_e) ~", v)), data = disc); b <- coef(f)[[1]]
  m0 <- mean(disc[[v]]); off <- b * (disc[[v]] - m0); tt <- disc$EP_ischemic_stroke_t; ee <- disc$EP_ischemic_stroke_e
  H0 <- sum(sapply(sort(unique(tt[ee == 1 & tt <= H])), function(u) sum(ee[tt == u]) / sum(exp(off[tt >= u]))))
  list(pred = 1 - exp(-H0 * exp(b * (rep[[v]] - m0))), slope_disc = b) }
PR <- list(SMART2 = risk3("SMART"), cfRNA = risk3("lock"), `SMART2 + cfRNA` = risk3("comb"))
for (k in names(PR)) rep[[paste0("p3_", k)]] <- PR[[k]]$pred
km_at <- function(idx) { f <- survfit(Surv(rep$EP_ischemic_stroke_t[idx], rep$EP_ischemic_stroke_e[idx]) ~ 1); s <- summary(f, times = H, extend = TRUE)
  c(obs = 1 - s$surv, lo = 1 - s$upper, hi = 1 - s$lower) }
CAL <- rbindlist(lapply(names(PR), function(k) { p <- rep[[paste0("p3_", k)]]; q <- cut(p, quantile(p, 0:3 / 3), include.lowest = TRUE, labels = FALSE)
  rbindlist(lapply(1:3, function(g) { o <- km_at(q == g); data.table(model = k, group = g, n = sum(q == g), pred = mean(p[q == g]), obs = o[["obs"]], lo = o[["lo"]], hi = o[["hi"]]) })) }))
ov <- km_at(rep(TRUE, nrow(rep)))
CS <- rbindlist(lapply(names(PR), function(k) { lp <- log(-log(1 - rep[[paste0("p3_", k)]]))
  cs <- coxph(Surv(EP_ischemic_stroke_t, EP_ischemic_stroke_e) ~ lp, data = cbind(rep, lp = lp))
  data.table(model = k, mean_pred = mean(rep[[paste0("p3_", k)]]), observed = ov[["obs"]], OE = ov[["obs"]] / mean(rep[[paste0("p3_", k)]]),
             cal_slope = coef(cs)[[1]], cal_slope_lo = confint(cs)[1], cal_slope_hi = confint(cs)[2]) }))
fwrite(CAL, file.path(OUT6, "fig6d_calibration_tertiles.csv")); fwrite(CS, file.path(OUT6, "fig6d_calibration_stats.csv"))
cat("[D] calibration at 3 y:\n"); print(CS[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
dd <- data.frame(t = rep$EP_ischemic_stroke_t, e = rep$EP_ischemic_stroke_e, SMART2 = rep$p3_SMART2, cfRNA = rep$p3_cfRNA, comb = rep$`p3_SMART2 + cfRNA`)
dc <- dca(Surv(t, e) ~ SMART2 + cfRNA + comb, data = dd, time = H, thresholds = seq(0, 0.15, by = 0.0025))
DC <- as.data.table(dc$dca)[, .(label = as.character(label), variable = as.character(variable), threshold, net_benefit)]
fwrite(DC, file.path(OUT6, "fig6e_dca.csv")); cat("[E] DCA rows:", nrow(DC), "\n")

## ---- F sex-specific HR per SD (both cohorts, both end points) -------------------------------------------
SXR <- list()
for (co in c("Discovery", "Replication")) { D <- if (co == "Discovery") copy(disc)[, score := lock] else copy(rep)[, score := lock]
  D[, z := (score - mean(score)) / sd(score)]
  for (ep in c("EP_major", "EP_ischemic_stroke")) { tt <- D[[paste0(ep, "_t")]]; ee <- D[[paste0(ep, "_e")]]
    pint <- summary(coxph(Surv(tt, ee) ~ z * male, data = D))$coefficients["z:male", 5]
    for (s in c(0, 1)) { k <- D$male == s; f <- coxph(Surv(tt[k], ee[k]) ~ z[k], data = D); ci <- exp(confint(f))
      rr <- pROC::roc(ee[k], D$score[k], direction = "<", quiet = TRUE)
      ac <- tryCatch(as.numeric(pROC::ci.auc(rr, method = "delong")), error = function(e) c(NA, NA, NA)); ac[2] <- as.numeric(pROC::auc(rr))
      if (any(!is.finite(ac[c(1, 3)])) || ac[1] == ac[3]) ac[c(1, 3)] <- NA
      SXR[[length(SXR) + 1]] <- data.table(cohort = co, endpoint = ep, sex = ifelse(s == 1, "Men", "Women"), n = sum(k), events = sum(ee[k]),
                                           HR = exp(coef(f))[[1]], lo = ci[1], hi = ci[2], p = summary(f)$coefficients[, 5], p_interaction = pint,
                                           AUC = ac[2], AUC_lo = ac[1], AUC_hi = ac[3]) } } }
SXR <- rbindlist(SXR); fwrite(SXR, file.path(OUT6, "fig6f_sex.csv"))
cat("[F] sex-specific HR per SD (events <10 suppressed):\n"); print(SXR[, .(cohort, endpoint, sex, n, events = sup(events), HR = round(HR, 2), AUC = round(AUC, 2), AUC_lo = round(AUC_lo, 2), AUC_hi = round(AUC_hi, 2), p_int = round(p_interaction, 3))])
cat("[ok] fig6 data written\n")

## ---- G incremental discrimination (IS, replication): C-index and 3-year AUC, SMART2 vs + cfRNA -------------
set.seed(20261002); NBS <- 1000
tt <- rep$EP_ischemic_stroke_t; ee <- rep$EP_ischemic_stroke_e; MV <- c(SMART2 = "SMART", cfRNA = "lock", `SMART2 + cfRNA` = "comb")
M0 <- cmask(tt, ee)
auc3 <- function(i, x) timeROC(T = tt[i], delta = ee[i], marker = x[i], cause = 1, weighting = "marginal", times = H, iid = FALSE)$AUC[[2]]
est <- sapply(MV, function(v) c(C = cidx(M0, rep[[v]]), AUC3 = auc3(seq_len(nrow(rep)), rep[[v]])))
BI <- replicate(NBS, sample.int(nrow(rep), replace = TRUE)); BI <- BI[, apply(BI, 2, function(j) sum(ee[j] * (tt[j] <= H)) >= 2)]
bt <- lapply(MV, function(v) t(apply(BI, 2, function(j) c(C = cidx(cmask(tt[j], ee[j]), rep[[v]][j]), AUC3 = auc3(j, rep[[v]])))))
INC <- rbindlist(lapply(names(MV), function(k) rbindlist(lapply(c("C", "AUC3"), function(m)
  data.table(model = k, metric = m, est = est[m, k], lo = quantile(bt[[k]][, m], .025, na.rm = TRUE), hi = quantile(bt[[k]][, m], .975, na.rm = TRUE))))))
d <- bt[["SMART2 + cfRNA"]] - bt[["SMART2"]]; pp <- function(v) 2 * min(mean(v <= 0, na.rm = TRUE), mean(v >= 0, na.rm = TRUE))
cmpA <- compare(timeROC(T = tt, delta = ee, marker = rep$SMART, cause = 1, weighting = "marginal", times = H, iid = TRUE),
                timeROC(T = tt, delta = ee, marker = rep$comb, cause = 1, weighting = "marginal", times = H, iid = TRUE))
DEL <- data.table(metric = c("C", "AUC3"), delta = c(est["C", "SMART2 + cfRNA"] - est["C", "SMART2"], est["AUC3", "SMART2 + cfRNA"] - est["AUC3", "SMART2"]),
                  lo = c(quantile(d[, "C"], .025, na.rm = TRUE), quantile(d[, "AUC3"], .025, na.rm = TRUE)),
                  hi = c(quantile(d[, "C"], .975, na.rm = TRUE), quantile(d[, "AUC3"], .975, na.rm = TRUE)),
                  p_boot = c(pp(d[, "C"]), pp(d[, "AUC3"])), p_timeROC = c(NA, as.numeric(cmpA$p_values_AUC[2])))
fwrite(INC, file.path(OUT6, "fig6g_increment.csv")); fwrite(DEL, file.path(OUT6, "fig6g_delta.csv"))
cat("[G] increment:\n"); print(INC[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); print(DEL[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

## ---- H decision-curve quantification at pre-specified thresholds ------------------------------------------------
## net benefit NB = TP/n - FP/n x pt/(1 - pt), TP/FP from the 3-year Kaplan-Meier risk (dcurves::dca, survival outcome)
## delta NB x 100 = additional net true positives per 100 patients vs the comparator;
## delta NB x (1 - pt)/pt x 100 = net reduction in unnecessary interventions per 100 patients at equal sensitivity
W6 <- dcast(DC[variable %in% c("SMART2", "comb", "all", "none")], threshold ~ variable, value.var = "net_benefit")
W6[, best_other := pmax(SMART2, all, none, na.rm = TRUE)]
Q <- W6[round(threshold, 4) %in% c(0.03, 0.05, 0.075, 0.10)][, .(threshold, NB_comb = comb, NB_SMART2 = SMART2, NB_all = all,
      dNB_vs_SMART2 = comb - SMART2, extraTP_per100_vs_SMART2 = 100 * (comb - SMART2),
      fewerFP_per100_vs_SMART2 = 100 * (comb - SMART2) * (1 - threshold) / threshold,
      dNB_vs_best_other = comb - best_other)]
W6 <- W6[order(threshold)][, better := !is.na(comb) & comb > best_other + 1e-6]
W6[, run := rleid(better)]
r5 <- W6[which.min(abs(threshold - 0.05)), run]                     # contiguous run containing 5%
win <- if (W6[run == r5, better[1]]) W6[run == r5, range(threshold)] else c(NA, NA)
fwrite(Q, file.path(OUT6, "fig6h_dca_quant.csv")); fwrite(data.table(from = win[1], to = win[2]), file.path(OUT6, "fig6h_dca_window.csv"))
cat(sprintf("[H] SMART2 + cfRNA has the highest net benefit for thresholds %.1f%% to %.1f%%\n", 100 * win[1], 100 * win[2]))
print(Q[, lapply(.SD, function(x) round(x, 4))])


## ---- I cumulative incidence (IS) by cohort: locked-model score > locked threshold (discovery top tertile) vs rest ---
disc[, group := grp(lock)]
KMd <- list(); RKd <- list(); KSd <- list()
for (co in c("Discovery", "Replication")) { D <- if (co == "Discovery") disc else rep
  tt <- D$EP_ischemic_stroke_t; ee <- D$EP_ischemic_stroke_e
  f <- survfit(Surv(tt, ee) ~ group, data = D); s1 <- summary(f, censored = TRUE)
  KMd[[co]] <- data.table(cohort = co, group = sub("group=", "", s1$strata), time = s1$time, cuminc = 1 - s1$surv, lo = 1 - s1$upper, hi = 1 - s1$lower)
  ar <- summary(f, times = 0:3, extend = TRUE); RKd[[co]] <- data.table(cohort = co, group = sub("group=", "", ar$strata), time = ar$time, n_risk = ar$n.risk)
  cx <- coxph(Surv(tt, ee) ~ group, data = D); lr <- survdiff(Surv(tt, ee) ~ group, data = D); ci <- exp(confint(cx))
  KSd[[co]] <- data.table(cohort = co, HR = exp(coef(cx))[[1]], lo = ci[1, 1], hi = ci[1, 2], p_logrank = 1 - pchisq(lr$chisq, 1)) }
fwrite(rbindlist(KMd), file.path(OUT6, "fig6i_km_IS_bycohort.csv")); fwrite(rbindlist(RKd), file.path(OUT6, "fig6i_atrisk_IS_bycohort.csv"))
fwrite(rbindlist(KSd), file.path(OUT6, "fig6i_stats_IS_bycohort.csv")); cat("[I] IS high vs low by cohort:\n"); print(rbindlist(KSd)[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])

## ---- J calibration curves at 3 years (IS, replication), three models: Cox model on cloglog(predicted risk), linear term
##      (df = 1; a spline was unstable for SMART2, whose predictions span a narrow range) (Austin et al., Stat Med 2020);
##      200-bootstrap 95% band; ICI = mean |observed - predicted|
library(splines); set.seed(20261003)
cal_curve <- function(p, t, e, grid) { x <- log(-log(1 - p)); xg <- log(-log(1 - grid))
  f <- coxph(Surv(t, e) ~ x, x = TRUE)
  obs_at <- function(xx) 1 - summary(survfit(f, newdata = data.frame(x = xx)), times = H, extend = TRUE)$surv
  list(grid = as.numeric(obs_at(xg)), ici = mean(abs(as.numeric(obs_at(x)) - p))) }
CC <- list(); ICI <- list()
for (k in c("SMART2", "cfRNA", "SMART2 + cfRNA")) { p <- rep[[paste0("p3_", k)]]; grid <- seq(quantile(p, .02), quantile(p, .98), length.out = 60)
  m <- cal_curve(p, rep$EP_ischemic_stroke_t, rep$EP_ischemic_stroke_e, grid)
  bs <- replicate(200, { j <- sample.int(length(p), replace = TRUE); tryCatch(cal_curve(p[j], rep$EP_ischemic_stroke_t[j], rep$EP_ischemic_stroke_e[j], grid)$grid, error = function(e) rep(NA, length(grid))) })
  CC[[k]] <- data.table(model = k, pred = grid, obs = m$grid, lo = apply(bs, 1, quantile, .025, na.rm = TRUE), hi = apply(bs, 1, quantile, .975, na.rm = TRUE))
  ICI[[k]] <- data.table(model = k, ICI = m$ici) }
PD <- rbindlist(lapply(c("SMART2", "cfRNA", "SMART2 + cfRNA"), function(k) data.table(model = k, pred = rep[[paste0("p3_", k)]])))
fwrite(rbindlist(CC), file.path(OUT6, "fig6j_calcurve.csv")); fwrite(rbindlist(ICI), file.path(OUT6, "fig6j_ici.csv"))
fwrite(PD[, .(n = .N), by = .(model, bin = round(pred / 0.005) * 0.005)], file.path(OUT6, "fig6j_pred_hist.csv"))
cat("[J] ICI:\n"); print(rbindlist(ICI))
