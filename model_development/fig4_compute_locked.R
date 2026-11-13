## ============================================================================
## Figure 4d data: locked five-gene model (refitted on the full discovery cohort,
## model_5gene.rds) evaluated on 8 end points in the discovery cohort (apparent).
##   risk score = frozen linear predictor; high risk = top tertile of discovery
##   scores (threshold saved for identical use in the replication cohort)
##   C-index from time to event; ROC-AUC, PR-AUC, sensitivity, specificity, NPV
##   from event occurrence during follow-up; 95% CI from 1,000 bootstrap resamples
## Output: 03_results/fig4/fig4d_metrics_locked.csv, locked_threshold.csv
## ============================================================================
suppressMessages({ library(data.table); library(haven); library(glmnet); library(survival)
                   library(pROC); library(PRROC); library(parallel) })
root <- "path/to/project/external"
OUT  <- file.path(root, "03_results", "fig4")
L    <- readRDS(file.path(root, "03_results", "disc344_to_external", "model_5gene.rds"))
EP <- c(`MACE` = "EP_major", `Ischaemic stroke` = "EP_ischemic_stroke", `Any stroke` = "EP_stroke",
        `Myocardial infarction` = "EP_MI", `Cardiovascular death` = "EP_CVdeath", `CAD events` = "EP_CAD",
        `Coronary events` = "EP_coronary", `Peripheral events` = "EP_peripheral")
D  <- fread("path/to/project/data/survival_data_cleanqc.csv", select = c("STUDY_NUMBER", L$sig))
cl <- read_sav(file.path(root, "01_raw_data", "data20240320.sav")); cl$STUDY_NUMBER <- as.integer(cl$STUDY_NUMBER)
C  <- cl[match(D$STUDY_NUMBER, cl$STUDY_NUMBER), ]
num <- function(x) { x <- as.numeric(zap_labels(x)); x[!is.na(x) & x < 0] <- NA; x }

lp  <- as.numeric(predict(L$fit, as.matrix(D[, L$sig, with = FALSE]), type = "link"))
thr <- unname(quantile(lp, 2 / 3)); high <- lp > thr
write.csv(data.frame(threshold_top_tertile = thr), file.path(OUT, "locked_threshold.csv"), row.names = FALSE)
cat("[ok] locked model applied | n =", length(lp), "| high-risk fraction =", round(mean(high), 3), "\n")

metrics <- function(t, e, x, h) {
  c(`C-index` = concordance(Surv(t, e) ~ x, reverse = TRUE)$concordance,
    `ROC-AUC` = as.numeric(auc(roc(e, x, direction = "<", quiet = TRUE))),
    `PR-AUC`  = pr.curve(scores.class0 = x[e == 1], scores.class1 = x[e == 0])$auc.integral,
    Sensitivity = mean(h[e == 1]), Specificity = mean(!h[e == 0]), NPV = sum(!h & e == 0) / sum(!h))
}
set.seed(20260929); BI <- replicate(1000, sample.int(length(lp), replace = TRUE))
res <- do.call(rbind, lapply(names(EP), function(k) {
  t <- num(C[[paste0(EP[[k]], "_time")]]); e <- num(C[[EP[[k]]]])
  ok <- !is.na(t) & !is.na(e); t <- t[ok]; e <- e[ok]; x <- lp[ok]; h <- high[ok]
  est <- metrics(t, e, x, h)
  bt <- do.call(rbind, mclapply(seq_len(ncol(BI)), function(b) {
    j <- BI[, b][BI[, b] <= length(x)]; if (sum(e[j]) < 2) return(rep(NA, 6)); metrics(t[j], e[j], x[j], h[j]) },
    mc.cores = max(1, detectCores() - 2)))
  data.frame(endpoint = k, n = length(e), events = sum(e), prevalence = mean(e), metric = names(est),
             estimate = unname(est), lo = apply(bt, 2, quantile, .025, na.rm = TRUE), hi = apply(bt, 2, quantile, .975, na.rm = TRUE))
}))
write.csv(res, file.path(OUT, "fig4d_metrics_locked.csv"), row.names = FALSE)
cat(sprintf("[check] MACE apparent C-index = %.3f (Table 2 apparent: 0.784)\n", res$estimate[res$endpoint == "MACE" & res$metric == "C-index"]))
cat("[ok] written: fig4d_metrics_locked.csv, locked_threshold.csv\n")
