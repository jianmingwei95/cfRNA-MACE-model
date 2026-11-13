## ============================================================================
## Strict external validation of "clinical score + cfRNA": combination weights locked in discovery
##   Combined LP = a x clinical score + b x cfRNA score. a, b = Cox coefficients estimated in the
##   discovery cohort (MACE, n = 344) and applied unchanged to the replication cohort (n = 277).
##   Clinical score internals = published coefficients (see benchmark_clinical_scores.R).
##   b estimated with the discovery out-of-sample (nested LOOCV) cfRNA score [primary] or with the
##   locked-model score (apparent) [sensitivity]; replication uses the locked-model score.
##   Replication metrics: C-index, AUC, PR-AUC (2,000-bootstrap CI); delta vs clinical score alone
##   (bootstrap P for C and PR-AUC, DeLong for AUC). End points: MACE and ischaemic stroke.
## Output: 03_results/benchmark_clinical/locked_weights_*.csv + Table_benchmark_locked_weights.xlsx
## ============================================================================
root <- "path/to/project/external"
src <- readLines(file.path(root, "clinical_benchmark", "benchmark_clinical_scores.R"))
eval(parse(text = src[seq_len(grep("^## ---------------- 0\\) cohorts", src) - 1)]))       # packages, metrics, scores(), impute()
suppressMessages({ library(glmnet) })
set.seed(42)

## ---- data --------------------------------------------------------------------------------------
lo <- fread("path/to/project/data/strict_nested_loocv_merged.csv", sep = ";", dec = ",")
Lk <- readRDS(file.path(root, "03_results", "disc344_to_external", "model_5gene.rds"))
Dg <- fread("path/to/project/data/survival_data_cleanqc.csv", select = c("STUDY_NUMBER", Lk$sig))
Dg[, lock := as.numeric(predict(Lk$fit, as.matrix(Dg[, Lk$sig, with = FALSE]), type = "link"))]
old <- fread("path/to/project/data/clinical_imp_final.csv")
disc <- merge(old[, .(STUDY_NUMBER)], lo[, .(STUDY_NUMBER, loocv = RiskScore)], by = "STUDY_NUMBER")
disc <- merge(disc, Dg[, .(STUDY_NUMBER, lock)], by = "STUDY_NUMBER")
rep <- fread(file.path(root, "03_results", "disc344_to_external", "external3batch_riskscore_and_genes.csv"))
SD <- scores(impute(disc$STUDY_NUMBER)$data); SR <- scores(impute(rep$STUDY_NUMBER)$data)
oD <- outcome(disc$STUDY_NUMBER, "EP_major")

## ---- 1) weights in discovery (MACE) --------------------------------------------------------------
W <- rbindlist(lapply(names(SD), function(k) rbindlist(lapply(c("loocv", "lock"), function(src) {
  f <- coxph(Surv(oD$t, oD$e) ~ clin + cf, data = data.frame(clin = SD[[k]], cf = disc[[src]]))
  data.table(score = k, cfRNA_for_b = src, a = coef(f)[["clin"]], b = coef(f)[["cf"]],
             a_se = sqrt(diag(vcov(f)))[1], b_se = sqrt(diag(vcov(f)))[2]) }))))
fwrite(W, file.path(OUT, "locked_weights_discovery.csv"))
cat("[weights] estimated in discovery (MACE):\n"); print(W[, .(score, cfRNA_for_b, a = round(a, 3), b = round(b, 3))])

## ---- 2) apply to replication -------------------------------------------------------------------
eval_fixed <- function(t, e, P, label) {                   # P: named list of fixed linear predictors
  n <- length(t); BI <- replicate(NB, sample.int(n, replace = TRUE)); BI <- BI[, apply(BI, 2, function(j) length(unique(e[j])) == 2)]
  M0 <- cmask(t, e)
  m3 <- function(j, x) { M <- if (is.null(j)) M0 else cmask(t[j], e[j]); jj <- if (is.null(j)) seq_len(n) else j
    c(C = cidx(M, x[jj]), AUC = aucf(e[jj], x[jj]), PR = apf(e[jj], x[jj])) }
  est <- lapply(P, function(x) m3(NULL, x)); bt <- lapply(P, function(x) t(apply(BI, 2, function(j) m3(j, x))))
  rbindlist(lapply(names(P), function(k) { ci <- apply(bt[[k]], 2, quantile, c(.025, .975), na.rm = TRUE)
    r <- data.table(analysis = label, model = k, n = n, events = sum(e),
      C = est[[k]][["C"]], C_lo = ci[1, "C"], C_hi = ci[2, "C"], AUC = est[[k]][["AUC"]], AUC_lo = ci[1, "AUC"], AUC_hi = ci[2, "AUC"],
      PR = est[[k]][["PR"]], PR_lo = ci[1, "PR"], PR_hi = ci[2, "PR"],
      dC = NA_real_, p_dC = NA_real_, dAUC = NA_real_, p_dAUC = NA_real_, dPR = NA_real_, p_dPR = NA_real_)
    if (grepl(" \\+ cfRNA$", k)) { b <- sub(" \\+ cfRNA$", "", k); d <- bt[[k]] - bt[[b]]
      pp <- function(v) 2 * min(mean(v <= 0, na.rm = TRUE), mean(v >= 0, na.rm = TRUE))
      r[, `:=`(dC = est[[k]][["C"]] - est[[b]][["C"]], p_dC = pp(d[, "C"]), dAUC = est[[k]][["AUC"]] - est[[b]][["AUC"]],
               p_dAUC = roc.test(roc(e, P[[b]], quiet = TRUE, direction = "<"), roc(e, P[[k]], quiet = TRUE, direction = "<"),
                                 method = "delong", paired = TRUE)$p.value,
               dPR = est[[k]][["PR"]] - est[[b]][["PR"]], p_dPR = pp(d[, "PR"]))] }
    r })) }
RES <- list()
for (src in c("loocv", "lock")) for (ep in c("EP_major", "EP_ischemic_stroke")) {
  o <- outcome(rep$STUDY_NUMBER, ep); P <- list()
  for (k in names(SR)) P[[k]] <- SR[[k]]
  P[["cfRNA"]] <- rep$risk
  for (k in names(SR)) { w <- W[score == k & cfRNA_for_b == src]; P[[paste0(k, " + cfRNA")]] <- w$a * SR[[k]] + w$b * rep$risk }
  lab <- sprintf("Replication | %s | weights locked in discovery (b from %s cfRNA score)", ifelse(ep == "EP_major", "MACE", "Ischaemic stroke"),
                 ifelse(src == "loocv", "LOOCV", "locked-model"))
  RES[[lab]] <- eval_fixed(o$t, o$e, P, lab); cat("[ok]", lab, "\n")
}
RES <- rbindlist(RES); fwrite(RES, file.path(OUT, "locked_weights_replication_all.csv"))

## ---- 3) publication-format tables -------------------------------------------------------------------
LAB <- c(SMART = "SMART", `SMART-REACH` = "SMART REACH", SMART2 = "SMART2", CEA = "CEA",
         `SMART + cfRNA` = "SMART + cfRNA", `SMART-REACH + cfRNA` = "SMART REACH + cfRNA", `SMART2 + cfRNA` = "SMART2 + cfRNA", `CEA + cfRNA` = "CEA + cfRNA")
ci  <- function(e, l, h) sprintf("%.2f (%.2f–%.2f)", e, l, h)
stars <- function(p) ifelse(p < 0.001, " ***", ifelse(p < 0.01, " **", ifelse(p < 0.05, " *", "")))
dlt <- function(d, p) ifelse(is.na(d), "–", paste0(sub("^-", "−", sprintf("%.2f", ifelse(abs(d) < 0.005, 0, d))), stars(p)))
wb <- createWorkbook(); modifyBaseFont(wb, fontSize = 10, fontName = "Arial")
hdr <- createStyle(textDecoration = "bold", border = "TopBottom", borderStyle = c("medium", "thin"), halign = "center", valign = "center", wrapText = TRUE)
FN <- c("Strict external validation: the Cox weights combining each clinical score with the cfRNA score (a x clinical + b x cfRNA) were estimated in the discovery cohort (MACE, n = 344) and applied unchanged to the replication cohort; no parameter was re-estimated in the replication cohort.",
        "Clinical scores use published coefficients; cfRNA = locked five-gene model. 95% CI from 2,000 bootstrap resamples.",
        "Δ: combined minus clinical score alone; ΔC-index and ΔPR-AUC by bootstrap, ΔAUC by DeLong. * P < 0.05, ** P < 0.01, *** P < 0.001.")
SH <- c(`Replication | MACE | weights locked in discovery (b from LOOCV cfRNA score)` = "Rep_MACE_locked",
        `Replication | Ischaemic stroke | weights locked in discovery (b from LOOCV cfRNA score)` = "Rep_IS_locked",
        `Replication | MACE | weights locked in discovery (b from locked-model cfRNA score)` = "Rep_MACE_locked_sens",
        `Replication | Ischaemic stroke | weights locked in discovery (b from locked-model cfRNA score)` = "Rep_IS_locked_sens")
for (a in names(SH)) { X <- RES[analysis == a][match(names(LAB), model)]; nm <- SH[[a]]
  T <- data.table(Model = LAB[X$model], `C-index (95% CI)` = ci(X$C, X$C_lo, X$C_hi), `AUC (95% CI)` = ci(X$AUC, X$AUC_lo, X$AUC_hi),
                  `PR-AUC (95% CI)` = ci(X$PR, X$PR_lo, X$PR_hi), `ΔC-index (p)` = dlt(X$dC, X$p_dC), `ΔAUC (p)` = dlt(X$dAUC, X$p_dAUC),
                  `ΔPR-AUC (p)` = dlt(X$dPR, X$p_dPR))
  addWorksheet(wb, nm, gridLines = FALSE); writeData(wb, nm, sprintf("%s (n = %d, %d events)", a, X$n[1], X$events[1]), startRow = 1)
  addStyle(wb, nm, createStyle(textDecoration = "bold", fontSize = 11), rows = 1, cols = 1)
  writeData(wb, nm, T, startRow = 3, headerStyle = hdr); addStyle(wb, nm, createStyle(halign = "center"), rows = 4:11, cols = 2:7, gridExpand = TRUE)
  addStyle(wb, nm, createStyle(border = "bottom", borderStyle = "thin"), rows = 7, cols = 1:7, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, nm, createStyle(border = "bottom", borderStyle = "medium"), rows = 11, cols = 1:7, gridExpand = TRUE, stack = TRUE)
  for (i in seq_along(FN)) { writeData(wb, nm, FN[i], startRow = 12 + i); mergeCells(wb, nm, cols = 1:7, rows = 12 + i) }
  addStyle(wb, nm, createStyle(fontSize = 8.5, wrapText = TRUE, valign = "top"), rows = 13:15, cols = 1, gridExpand = TRUE)
  setRowHeights(wb, nm, rows = 13:15, heights = 28); setColWidths(wb, nm, 1:7, c(22, 18, 18, 18, 14, 13, 14))
  cat("\n##", a, "\n"); print(T, row.names = FALSE) }
addWorksheet(wb, "Weights_discovery"); writeData(wb, "Weights_discovery", W)
saveWorkbook(wb, file.path(OUT, "Table_benchmark_locked_weights.xlsx"), overwrite = TRUE)
cat("\n[ok] written Table_benchmark_locked_weights.xlsx\n")
