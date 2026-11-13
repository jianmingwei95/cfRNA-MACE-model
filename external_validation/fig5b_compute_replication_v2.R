## ============================================================================
## Figure 4f data: locked five-gene model in the replication cohort (serum, n = 277)
##   risk score = frozen linear predictor (external3batch_riskscore_and_genes.csv)
##   high risk = discovery top-tertile threshold (locked_threshold.csv), no re-estimation
##   same metrics and 1,000-bootstrap CIs as fig4_compute_locked.R (discovery)
## v2 (2026-10-02): adds PPV. Output: 03_results/fig5/fig5b_metrics_replication_v2.csv (aggregate only)
## ============================================================================
suppressMessages({ library(data.table); library(haven); library(survival); library(pROC); library(PRROC); library(parallel) })
root <- "path/to/project/external"
OUT  <- file.path(root, "03_results", "fig5"); dir.create(OUT, showWarnings = FALSE)
IN4  <- file.path(root, "03_results", "fig4")
EP <- c(`MACE` = "EP_major", `Ischaemic stroke` = "EP_ischemic_stroke", `Any stroke` = "EP_stroke",
        `Myocardial infarction` = "EP_MI", `Cardiovascular death` = "EP_CVdeath", `CAD events` = "EP_CAD",
        `Coronary events` = "EP_coronary", `Peripheral events` = "EP_peripheral")
sc <- fread(file.path(root, "03_results", "disc344_to_external", "external3batch_riskscore_and_genes.csv"))
cl <- read_sav(file.path(root, "01_raw_data", "data20240320.sav")); cl$STUDY_NUMBER <- as.integer(cl$STUDY_NUMBER)
C  <- cl[match(sc$STUDY_NUMBER, cl$STUDY_NUMBER), ]
num <- function(x) { x <- as.numeric(zap_labels(x)); x[!is.na(x) & x < 0] <- NA; x }
thr <- read.csv(file.path(IN4, "locked_threshold.csv"))$threshold_top_tertile
lp <- sc$risk; high <- lp > thr
stopifnot(nrow(sc) == 277, !anyNA(lp))
cat(sprintf("[ok] n = %d | discovery threshold = %.3f | high-risk fraction in replication = %.3f\n", length(lp), thr, mean(high)))
cat(sprintf("[scale] replication risk score quantiles 33/67%%: %.3f / %.3f\n", quantile(lp, 1/3), quantile(lp, 2/3)))

metrics <- function(t, e, x, h) {
  c(`C-index` = concordance(Surv(t, e) ~ x, reverse = TRUE)$concordance,
    `ROC-AUC` = as.numeric(auc(roc(e, x, direction = "<", quiet = TRUE))),
    `PR-AUC`  = pr.curve(scores.class0 = x[e == 1], scores.class1 = x[e == 0])$auc.integral,
    Sensitivity = mean(h[e == 1]), Specificity = mean(!h[e == 0]), PPV = sum(h & e == 1) / sum(h), NPV = sum(!h & e == 0) / sum(!h))
}
set.seed(20260929); BI <- replicate(1000, sample.int(length(lp), replace = TRUE))
res <- do.call(rbind, lapply(names(EP), function(k) {
  t <- num(C[[paste0(EP[[k]], "_time")]]); e <- num(C[[EP[[k]]]])
  ok <- !is.na(t) & !is.na(e); t <- t[ok]; e <- e[ok]; x <- lp[ok]; h <- high[ok]
  est <- metrics(t, e, x, h)
  bt <- do.call(rbind, mclapply(seq_len(ncol(BI)), function(b) {
    j <- BI[, b][BI[, b] <= length(x)]; if (sum(e[j]) < 2) return(rep(NA, 7)); metrics(t[j], e[j], x[j], h[j]) },
    mc.cores = 3))
  data.frame(endpoint = k, n = length(e), events = sum(e), prevalence = mean(e), metric = names(est),
             estimate = unname(est), lo = apply(bt, 2, quantile, .025, na.rm = TRUE), hi = apply(bt, 2, quantile, .975, na.rm = TRUE))
}))
write.csv(res, file.path(OUT, "fig5b_metrics_replication_v2.csv"), row.names = FALSE)
print(data.table(res)[metric == "C-index", .(endpoint, events, C = round(estimate, 3))])
cat("[ok] written: fig5b_metrics_replication_v2.csv\n")
