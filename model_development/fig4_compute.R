## ============================================================================
## Figure 4 data: model development and finalisation (discovery cohort, n = 344)
##   a  LOOCV selection frequency of the top-10 genes
##   b  univariable Cox HR per SD for MACE, top-10 genes
##   d  fixed five-gene model, LOOCV (lambda re-tuned in every training set),
##      MACE-trained cross-validated risk score evaluated on 8 end points:
##      C-index (time-to-event), ROC-AUC and PR-AUC (event during follow-up),
##      sensitivity / specificity / NPV at an in-fold threshold (training-set
##      67th percentile = top tertile), 95% CI by 1,000 bootstrap resamples
## Outputs (aggregate only): 03_results/fig4/*.csv ; CV scores kept in .rds
## ============================================================================
suppressMessages({ library(data.table); library(haven); library(glmnet); library(survival)
                   library(pROC); library(PRROC); library(parallel) })
root <- "path/to/project/external"
OUT  <- file.path(root, "03_results", "fig4"); dir.create(OUT, showWarnings = FALSE)
TOP10 <- c("EI24", "EIF4E", "SPINT2", "TRIM25", "F11R", "RFXANK", "NR2C1", "ATP1A1", "TCAF1", "TSC22D4")
SIG   <- c("SPINT2", "EIF4E", "TRIM25", "EI24", "F11R")
EP <- c(`MACE` = "EP_major", `Ischaemic stroke` = "EP_ischemic_stroke", `Any stroke` = "EP_stroke",
        `Myocardial infarction` = "EP_MI", `Cardiovascular death` = "EP_CVdeath", `CAD events` = "EP_CAD",
        `Coronary events` = "EP_coronary", `Peripheral events` = "EP_peripheral")
CORES <- max(1, detectCores() - 1)

D  <- fread("path/to/project/data/survival_data_cleanqc.csv", select = c("STUDY_NUMBER", "EP_major", "EP_major_time", TOP10))
cl <- read_sav(file.path(root, "01_raw_data", "data20240320.sav")); cl$STUDY_NUMBER <- as.integer(cl$STUDY_NUMBER)
C  <- cl[match(D$STUDY_NUMBER, cl$STUDY_NUMBER), ]
num <- function(x) { x <- as.numeric(zap_labels(x)); x[!is.na(x) & x < 0] <- NA; x }
stopifnot(nrow(D) == 344, all(num(C$EP_major) == D$EP_major, na.rm = TRUE))
cat("[ok] discovery n =", nrow(D), "\n")

## ---- a ----------------------------------------------------------------------
st <- read.csv("path/to/project/data/strict_nested_loocv_10features_feature_stability.csv")
st <- st[order(-st$Count), ][1:10, ]; stopifnot(identical(st$Feature, TOP10))
write.csv(data.frame(gene = st$Feature, count = st$Count, pct = 100 * st$Count / 344, locked = st$Feature %in% SIG),
          file.path(OUT, "fig4a_selection.csv"), row.names = FALSE)

## ---- b ----------------------------------------------------------------------
hb <- do.call(rbind, lapply(TOP10, function(g) {
  s <- summary(coxph(Surv(D$EP_major_time, D$EP_major) ~ scale(D[[g]])[, 1]))
  data.frame(gene = g, HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5])
}))
hb$locked <- hb$gene %in% SIG
write.csv(hb, file.path(OUT, "fig4b_univariable_hr.csv"), row.names = FALSE)
cat("[ok] b: top-5 all HR > 1:", all(hb$HR[hb$locked] > 1), "\n")

## ---- d: LOOCV with fixed five genes ----------------------------------------------
X <- as.matrix(D[, ..SIG]); tt <- D$EP_major_time; ee <- D$EP_major
cvres <- mclapply(seq_len(nrow(X)), function(i) {
  set.seed(1000 + i)
  cv  <- cv.glmnet(X[-i, ], Surv(tt[-i], ee[-i]), family = "cox", alpha = 0, nfolds = 10)
  lpt <- as.numeric(predict(cv, X[-i, ], s = "lambda.min"))
  c(lp = as.numeric(predict(cv, X[i, , drop = FALSE], s = "lambda.min")), thr = unname(quantile(lpt, 2 / 3)))
}, mc.cores = CORES)
cvres <- do.call(rbind, cvres); lp <- cvres[, "lp"]; high <- lp > cvres[, "thr"]
saveRDS(list(STUDY_NUMBER = D$STUDY_NUMBER, lp = lp, high = high), file.path(OUT, "fig4_cv_scores.rds"))
cat("[ok] LOOCV done | high-risk fraction =", round(mean(high), 3), "\n")

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
    mc.cores = CORES))
  data.frame(endpoint = k, n = length(e), events = sum(e), prevalence = mean(e), metric = names(est),
             estimate = unname(est), lo = apply(bt, 2, quantile, .025, na.rm = TRUE), hi = apply(bt, 2, quantile, .975, na.rm = TRUE))
}))
write.csv(res, file.path(OUT, "fig4d_metrics.csv"), row.names = FALSE)
cat("[check] events per end point:", paste(unique(res[, c("endpoint", "events")])$events, collapse = ","), "\n")
cat(sprintf("[check] MACE CV C-index = %.3f (previous 10-fold CV estimate 0.672)\n", res$estimate[res$endpoint == "MACE" & res$metric == "C-index"]))
cat("[ok] written to 03_results/fig4\n")
