## ============================================================================
## Clinical risk scores (SMART, SMART-REACH, SMART2, CEA) vs + cfRNA five-gene score
## Applied to discovery (n = 344) and replication (n = 277).
##
## Evaluation: Cox model on each clinical LP, and on LP + cfRNA score, fitted
## in-sample; predicted partial hazard -> Harrell C (time-to-event), ROC-AUC and PR-AUC (average
## precision) on the event indicator; 95% CI from 2,000 bootstrap resamples of patients;
## delta C and delta PR-AUC: bootstrap, two-sided P = 2 x min(P(d <= 0), P(d >= 0)); delta AUC: DeLong.
##
## Clinical scores use the published coefficients. The 'years since first ASCVD diagnosis' term of
## SMART and SMART2 is omitted because the diagnosis date is not available in this cohort.
## hsCRP <= 0 or > 300 mg/L treated as implausible (set to missing, then imputed); applies to both cohorts.
## cfRNA score: discovery = strict nested LOOCV RiskScore;
##              replication = locked five-gene model (coefficients fixed).
## End points: MACE (EP_major) and ischaemic stroke.
## Outputs (aggregate only): 03_results/benchmark_clinical/*.csv + Table_benchmark_clinical.xlsx
## ============================================================================
suppressMessages({ library(haven); library(mice); library(survival); library(pROC); library(data.table); library(openxlsx) })
root <- "path/to/project/external"
OUT <- file.path(root, "03_results", "benchmark_clinical"); dir.create(OUT, showWarnings = FALSE)
NB <- 2000; set.seed(42)
cl <- read_sav(file.path(root, "01_raw_data", "data20240320.sav")); cl$STUDY_NUMBER <- as.integer(cl$STUDY_NUMBER)
num <- function(x) { x <- as.numeric(zap_labels(x)); x[!is.na(x) & x %in% c(-999, -888)] <- NA; x }

## ---------------- metrics ------------------------------------------------------------------------
cmask <- function(t, e) (e == 1) & (outer(t, t, "<") | (outer(t, t, "==") & outer(e == 1, e == 0, "&")))
cidx  <- function(M, x) { d <- outer(x, x, "-"); (sum(M & d > 0) + 0.5 * sum(M & d == 0)) / sum(M) }
aucf  <- function(y, x) { r <- rank(x); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
apf   <- function(y, x) { o <- order(-x); ys <- y[o]; xs <- x[o]; last <- c(xs[-1] != xs[-length(xs)], TRUE)
  tp <- cumsum(ys)[last]; fp <- cumsum(1 - ys)[last]; prec <- tp / (tp + fp); rec <- tp / sum(y); sum(diff(c(0, rec)) * prec) }
coxlp <- function(t, e, X) { f <- coxph(Surv(t, e) ~ ., data = as.data.frame(X)); as.numeric(predict(f, type = "lp")) }
evaluate <- function(t, e, clin, risk, label) {
  n <- length(t); BI <- replicate(NB, sample.int(n, replace = TRUE)); ok <- apply(BI, 2, function(j) length(unique(e[j])) == 2)
  BI <- BI[, ok]; M0 <- cmask(t, e)
  m3 <- function(j, x) { M <- if (is.null(j)) M0 else cmask(t[j], e[j]); jj <- if (is.null(j)) seq_len(n) else j
    c(C = cidx(M, x[jj]), AUC = aucf(e[jj], x[jj]), PR = apf(e[jj], x[jj])) }
  preds <- list(); for (k in names(clin)) preds[[k]] <- coxlp(t, e, data.frame(lp = clin[[k]]))
  preds[["cfRNA"]] <- coxlp(t, e, data.frame(risk = risk))
  for (k in names(clin)) preds[[paste0(k, " + cfRNA")]] <- coxlp(t, e, data.frame(lp = clin[[k]], risk = risk))
  est <- lapply(preds, function(x) m3(NULL, x))
  bt  <- lapply(preds, function(x) t(apply(BI, 2, function(j) m3(j, x))))
  rows <- rbindlist(lapply(names(preds), function(k) {
    ci <- apply(bt[[k]], 2, quantile, c(.025, .975), na.rm = TRUE)
    r <- data.table(analysis = label, model = k, n = n, events = sum(e),
                    C = est[[k]][["C"]], C_lo = ci[1, "C"], C_hi = ci[2, "C"],
                    AUC = est[[k]][["AUC"]], AUC_lo = ci[1, "AUC"], AUC_hi = ci[2, "AUC"],
                    PR = est[[k]][["PR"]], PR_lo = ci[1, "PR"], PR_hi = ci[2, "PR"],
                    dC = NA_real_, p_dC = NA_real_, dAUC = NA_real_, p_dAUC = NA_real_, dPR = NA_real_, p_dPR = NA_real_)
    if (grepl(" \\+ cfRNA$", k)) { b <- sub(" \\+ cfRNA$", "", k)
      d <- bt[[k]] - bt[[b]]; pp <- function(v) 2 * min(mean(v <= 0, na.rm = TRUE), mean(v >= 0, na.rm = TRUE))
      r[, `:=`(dC = est[[k]][["C"]] - est[[b]][["C"]], p_dC = pp(d[, "C"]),
               dAUC = est[[k]][["AUC"]] - est[[b]][["AUC"]],
               p_dAUC = roc.test(roc(e, preds[[b]], quiet = TRUE, direction = "<"), roc(e, preds[[k]], quiet = TRUE, direction = "<"),
                                 method = "delong", paired = TRUE)$p.value,
               dPR = est[[k]][["PR"]] - est[[b]][["PR"]], p_dPR = pp(d[, "PR"]))] }
    r }))
  rows
}

## ---------------- score construction (published coefficients) --------------------------------------
scores <- function(D) {
  with(D, {
    Male <- as.numeric(Gender == 1); Smoking <- SmokerCurrent; Diabetes <- DM
    NonHDL <- TC_all - HDL_all
    SMART2 <- 0.330356631 * heart325 - 0.03496022 * Age + 0.000551072 * Age^2 - 0.21072103 * Med.acetylsal +
      0.294701954 * CAD_history - 0.03967521 * GFR_MDRD + 0.000218613 * GFR_MDRD^2 + 0.151760173 * log(hsCRP_plasma) +
      0.345583271 * Smoking + 0.34831786 * StrokeTIA_history + 0.318170659 * Diabetes + 0.540364249 * log(NonHDL) +
      0.22446658 * Peripheral.interv + 0.018913154 * systolic + 0.287658743 * Male
    vl <- CAD_history + StrokeTIA_history + Peripheral.interv + heart325
    REACH <- 0.0720 * Male + 0.4309 * Smoking + 0.4357 * Diabetes - 0.0281 * systolic + 0.0001 * systolic^2 -
      0.3671 * TC_all + 0.0356 * TC_all^2 + 0.0061 * creat + 0.3176 * (vl == 2) + 0.2896 * (vl >= 3)
    SMART <- -0.0850 * Age + 0.0011 * Age^2 + 0.1561 * Male + 0.2232 * Diabetes + 0.2617 * Smoking + 0.0043 * (systolic / 10) +
      0.0959 * TC_all - 0.4256 * HDL_all + 0.1394 * log(hsCRP_plasma) - 0.0532 * GFR_MDRD + 0.0003 * GFR_MDRD^2 +
      0.4058 * StrokeTIA_history + 0.1401 * CAD_history + 0.5578 * heart325 + 0.2832 * Peripheral.interv
    I_Either <- as.numeric(CADPAOD_history == 1); I_Both <- as.numeric(CADPAOD_history == 2)
    I_TIA <- as.numeric(Symptoms.4g %in% c(1, 2)); I_Stroke <- as.numeric(Symptoms.4g == 3); Contra <- stenosis_con_bin
    CEA <- 0.0405 * Age + 0.0734 * I_Either + 0.4387 * I_Both + 0.2747 * Smoking - 0.0050 * systolic +
      0.3685 * Hypertension.drugs + 0.3520 * I_TIA + 0.5026 * I_Stroke + 0.3338 * Contra + 0.5035 * log(creat)
    list(SMART = SMART, `SMART-REACH` = REACH, SMART2 = SMART2, CEA = CEA)
  })
}

## ---------------- imputation (MICE: pmm / logreg / polyreg, seed 2026, m = 5, 1st set) --
VB <- c("heart325", "Med.acetylsal", "CAD_history", "SmokerCurrent", "StrokeTIA_history", "DM", "Peripheral.interv",
        "Hypertension.drugs", "stenosis_con_bin")
VC <- c("GFR_MDRD", "hsCRP_plasma", "TC_all", "HDL_all", "systolic", "creat")
impute <- function(ids) {
  E <- cl[match(ids, cl$STUDY_NUMBER), ]
  D <- data.frame(Age = num(E$Age), Gender = num(E$Gender), EP_major = num(E$EP_major), EP_major_time = num(E$EP_major_time))
  for (v in c(VB, VC)) D[[v]] <- num(E[[v]])
  D$hsCRP_plasma[!is.na(D$hsCRP_plasma) & (D$hsCRP_plasma <= 0 | D$hsCRP_plasma > 300)] <- NA   # implausible -> impute
  D$Symptoms.4g <- factor(num(E$Symptoms.4g), levels = 0:3)
  D$CADPAOD_history <- factor(num(E$CADPAOD_history), levels = 0:2)
  D$CADPAOD_2level <- factor(num(E$CADPAOD_history), levels = c(0, 2))               # auxiliary two-level coding; not a predictor, not used in the scores
  for (v in VB) D[[v]] <- factor(D[[v]], levels = 0:1)
  meth <- make.method(D); meth[VC] <- "pmm"; meth[VB] <- "logreg"; meth[c("Symptoms.4g", "CADPAOD_history")] <- "polyreg"
  meth["CADPAOD_2level"] <- "logreg"; meth[c("Age", "Gender", "EP_major", "EP_major_time")] <- ""
  pred <- make.predictorMatrix(D); diag(pred) <- 0; pred[, "CADPAOD_2level"] <- 0; pred["CADPAOD_2level", "CADPAOD_history"] <- 0
  set.seed(2026); imp <- mice(D, m = 5, method = meth, predictorMatrix = pred, maxit = 20, printFlag = FALSE)
  C <- complete(imp, 1)
  for (v in c(VB, "Symptoms.4g", "CADPAOD_history", "CADPAOD_2level")) C[[v]] <- as.numeric(as.character(C[[v]]))
  list(data = C, n_missing = colSums(is.na(D)))
}
outcome <- function(ids, ep) { E <- cl[match(ids, cl$STUDY_NUMBER), ]
  list(t = num(E[[paste0(ep, "_time")]]), e = as.numeric(num(E[[ep]]) %in% 1)) }

## ---------------- 0) cohorts ----------------
old <- fread("path/to/project/data/clinical_imp_final.csv"); lo <- fread("path/to/project/data/strict_nested_loocv_merged.csv", sep = ";", dec = ",")
if (!is.numeric(lo$RiskScore)) lo[, RiskScore := as.numeric(gsub(",", ".", RiskScore))]
OD <- merge(old[, .(STUDY_NUMBER)], lo[, .(STUDY_NUMBER, RiskScore)], by = "STUDY_NUMBER")

## ---------------- 1) impute + score both cohorts, two end points ------------------------
disc_ids <- OD$STUDY_NUMBER
rep <- fread(file.path(root, "03_results", "disc344_to_external", "external3batch_riskscore_and_genes.csv"))
COH <- list(Discovery = list(ids = disc_ids, risk = OD$RiskScore), Replication = list(ids = rep$STUDY_NUMBER, risk = rep$risk))
ALL <- list(); MISS <- list()
for (co in names(COH)) {
  im <- impute(COH[[co]]$ids); MISS[[co]] <- im$n_missing
  S <- scores(im$data)
  for (ep in c("EP_major", "EP_ischemic_stroke")) {
    o <- outcome(COH[[co]]$ids, ep)
    lab <- sprintf("%s | %s", co, ifelse(ep == "EP_major", "MACE", "Ischaemic stroke"))
    ALL[[lab]] <- evaluate(o$t, o$e, S, COH[[co]]$risk, lab); cat("[ok]", lab, "\n")
  }
}
RES <- rbindlist(ALL); fwrite(RES, file.path(OUT, "benchmark_clinical_all.csv"))
fwrite(as.data.table(do.call(cbind, MISS), keep.rownames = "variable"), file.path(OUT, "benchmark_missing_before_imputation.csv"))

## ---------------- 2) formatted tables (one sheet per analysis) -----------------------------------------------
f2 <- function(e, l, h) sprintf("%.2f (%.2f–%.2f)", e, l, h)
star <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, " ***", ifelse(p < 0.01, " **", ifelse(p < 0.05, " *", ""))))
fd <- function(d, p) ifelse(is.na(d), "–", sprintf("%+.2f%s", d, star(p)))
wb <- createWorkbook(); hdr <- createStyle(textDecoration = "bold", border = "TopBottom", borderStyle = "medium", wrapText = TRUE, valign = "center")
for (a in unique(RES$analysis)) {
  X <- RES[analysis == a]; ord <- c("SMART", "SMART-REACH", "SMART2", "CEA", "cfRNA", paste(c("SMART", "SMART-REACH", "SMART2", "CEA"), "+ cfRNA"))
  X <- X[match(ord, model)]
  T <- data.table(Model = X$model, `C-index (95% CI)` = f2(X$C, X$C_lo, X$C_hi), `AUC (95% CI)` = f2(X$AUC, X$AUC_lo, X$AUC_hi),
                  `PR-AUC (95% CI)` = f2(X$PR, X$PR_lo, X$PR_hi), `ΔC-index` = fd(X$dC, X$p_dC), `ΔAUC` = fd(X$dAUC, X$p_dAUC),
                  `ΔPR-AUC` = fd(X$dPR, X$p_dPR))
  p <- trimws(strsplit(a, "\\|")[[1]]); sh <- paste(substr(p[1], 1, 4), ifelse(p[2] == "MACE", "MACE", "IS"), sep = "_")
  addWorksheet(wb, sh)
  writeData(wb, sh, sprintf("%s (n = %d, %d events)", a, X$n[1], X$events[1]), startRow = 1)
  writeData(wb, sh, T, startRow = 3, headerStyle = hdr); setColWidths(wb, sh, 1:7, c(22, 20, 20, 20, 12, 12, 12))
  writeData(wb, sh, "* P < 0.05, ** P < 0.01, *** P < 0.001 (delta C and delta PR-AUC: bootstrap; delta AUC: DeLong). All Cox models fitted in-sample.", startRow = 14)
}
saveWorkbook(wb, file.path(OUT, "Table_benchmark_clinical.xlsx"), overwrite = TRUE)
cat("[ok] written benchmark_clinical_all.csv, Table_benchmark_clinical.xlsx\n")
