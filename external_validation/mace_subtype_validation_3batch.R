## ============================================================================
## MACE-subtype validation of the 5-gene cfRNA signature
##   SPINT2, EIF4E, TRIM25, EI24, F11R
##
## Discovery (train, n = 344) -> independent external, 3 BATCHES
##   MOK12069 + MOK12755 + MOK12754
## Frozen ridge-Cox trained on MACE (EP_major), applied WITHOUT refitting.
##
## Pipeline: samples that are also in the discovery cohort are removed first; then sample QC,
## gene QC, quantile normalisation, ComBat (batch only), reference quantile mapping onto the
## discovery distribution, and scoring with the locked model. Eight end points are evaluated.
## ============================================================================
suppressMessages({
  library(org.Hs.eg.db); library(AnnotationDbi); library(limma); library(sva)
  library(haven); library(dplyr); library(glmnet); library(survival); library(survivalROC)
})

root <- "path/to/project/external"
raw  <- file.path(root, "01_raw_data")
resd <- file.path(root, "03_results", "disc344_to_external")
SIG  <- c("SPINT2", "EIF4E", "TRIM25", "EI24", "F11R")

ENDPT <- list(
  `Ischemic stroke`       = c("EP_ischemic_stroke", "EP_ischemic_stroke_time"),
  `Stroke (any)`          = c("EP_stroke",          "EP_stroke_time"),
  `Myocardial infarction` = c("EP_MI",              "EP_MI_time"),
  `CV death`              = c("EP_CVdeath",         "EP_CVdeath_time"),
  `MACE (composite)`      = c("EP_major",           "EP_major_time"),
  `CAD`                   = c("EP_CAD",             "EP_CAD_time"),
  `Peripheral`            = c("EP_peripheral",      "EP_peripheral_time"),
  `Coronary`              = c("EP_coronary",        "EP_coronary_time"))

Cidx <- function(t, s, r) {
  ok <- !is.na(t) & !is.na(s) & !is.na(r)
  if (sum(s[ok] %in% 1) < 3) return(c(NA, NA))
  summary(coxph(Surv(t[ok], s[ok]) ~ r[ok]))$concordance
}
auc3 <- function(t, s, r) {
  ok <- !is.na(t) & !is.na(s) & !is.na(r)
  if (sum(s[ok] %in% 1) < 3) return(NA)
  tryCatch(survivalROC(t[ok], s[ok], r[ok], predict.time = 3, method = "KM")$AUC,
           error = function(e) NA)
}
hr_sd <- function(t, s, r) {
  ok <- !is.na(t) & !is.na(s) & !is.na(r)
  if (sum(s[ok] %in% 1) < 3) return(c(NA, NA, NA, NA))
  cf <- summary(coxph(Surv(t[ok], s[ok]) ~ scale(r[ok])[, 1]))
  c(cf$conf.int[1, 1], cf$conf.int[1, 3], cf$conf.int[1, 4], cf$coefficients[1, 5])
}

## ---- 0. read raw counts of the three batches ---------------------------------------------------
c1 <- read.table(file.path(raw, "cfRNAs_MOK12069_combined_raw_counts.txt"), header = TRUE, row.names = 1, check.names = FALSE)
c2 <- read.table(file.path(raw, "cfRNAs_MOK12755_combined_raw_counts.txt"), header = TRUE, row.names = 1, check.names = FALSE)
c3 <- read.table(file.path(raw, "cfRNAs_MOK12754_combined_raw_counts.txt"), header = TRUE, row.names = 1, check.names = FALSE)
common <- Reduce(intersect, list(rownames(c1), rownames(c2), rownames(c3)))
CNT <- cbind(c1[common, , drop = FALSE], c2[common, , drop = FALSE], c3[common, , drop = FALSE])
BAT <- factor(c(rep("MOK12069", ncol(c1)), rep("MOK12755", ncol(c2)), rep("MOK12754", ncol(c3))),
              levels = c("MOK12069", "MOK12755", "MOK12754"))
CNT <- CNT[grep("^ERCC", rownames(CNT), invert = TRUE), ]
colnames(CNT) <- gsub("\\D+", "", colnames(CNT)); names(BAT) <- colnames(CNT)
stopifnot(!anyDuplicated(colnames(CNT)))
cat("3-batch raw (ERCC removed):", dim(CNT), "| A", ncol(c1), "B", ncol(c2), "C", ncol(c3), "\n")

## ---- 1. remove samples that are also in the discovery cohort before any QC ------
disc <- read.csv("path/to/project/data/survival_data_cleanqc.csv", check.names = FALSE)
disc_id <- as.integer(disc$STUDY_NUMBER)
ov <- intersect(as.integer(colnames(CNT)), disc_id)
cat("samples overlapping with discovery (removed before QC):", length(ov), "\n")
CNT <- CNT[, !(as.integer(colnames(CNT)) %in% ov), drop = FALSE]
BAT <- droplevels(BAT[colnames(CNT)])
cat("replication samples entering QC:", ncol(CNT), "|",
    paste(names(table(BAT)), table(BAT), collapse = " "), "\n")

clin <- read_sav(file.path(raw, "data20240320.sav")); clin$STUDY_NUMBER <- as.integer(clin$STUDY_NUMBER)
L <- readRDS(file.path(resd, "model_5gene.rds")); stopifnot(identical(L$sig, SIG))
cat("\n== frozen coefficients (NOT refitted) ==\n")
print(round(as.matrix(coef(L$fit))[SIG, , drop = FALSE], 4))

qmap <- function(x, ref) {
  p <- (rank(x, ties.method = "average") - 0.5) / length(x)
  as.numeric(quantile(ref, probs = p, type = 7, names = FALSE))
}

## ---- 2. QC -> normalisation -> ComBat -> locked score -------------------------
##   sample QC: library size >= 5th percentile and > 10,000 detected genes
##   gene QC  : CPM > 0.5 in >= 70% of samples
GD_THR   <- 10000
GENE_FRC <- 0.70

build <- function(frac, tag) {
  cnt <- CNT; bat <- BAT
  tot <- colSums(cnt); k <- tot >= quantile(tot, 0.05) & colSums(cnt > 0) > GD_THR
  cnt <- cnt[, k]; bat <- droplevels(bat[k])
  libs <- colSums(cnt); cpm <- sweep(cnt, 2, libs, "/") * 1e6
  cnt <- cnt[rowSums(cpm > 0.5) >= frac * ncol(cnt), ]
  sym <- mapIds(org.Hs.eg.db, rownames(cnt), "SYMBOL",   "ENSEMBL", multiVals = "first")
  gt  <- mapIds(org.Hs.eg.db, rownames(cnt), "GENETYPE", "ENSEMBL", multiVals = "first")
  ok  <- !is.na(sym) & gt %in% "protein-coding" & !duplicated(sym)
  cnt <- cnt[ok, ]; rownames(cnt) <- sym[ok]
  cnt <- cnt[!grepl("^RP[LS]", rownames(cnt)), ]
  missing <- setdiff(SIG, rownames(cnt))
  cat(sprintf("\n[%s] after sample QC (detected genes > %d) n = %d | genes after filter (%.0f%%) = %d | missing signature genes: %s\n",
              tag, GD_THR, ncol(cnt), 100*frac, nrow(cnt),
              if (length(missing)) paste(missing, collapse=",") else "none"))
  if (length(missing)) {
    cat("   -> signature incomplete, cannot score\n")
    return(NULL)
  }
  libs <- colSums(cnt); cpm <- sweep(cnt, 2, libs, "/") * 1e6
  m <- as.matrix(log2(normalizeQuantiles(cpm) + 1))
  m <- m[apply(m, 1, function(g) all(tapply(g, bat, function(v) sd(v) > 0))), ]
  ext <- ComBat(dat = m, batch = bat, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
  X <- sapply(SIG, function(g) qmap(ext[g, ], disc[[g]]))
  list(id = as.integer(colnames(ext)), bat = bat, ngene = nrow(ext),
       risk = as.numeric(predict(L$fit, X[, SIG, drop = FALSE], type = "link")), X = X)
}

CFG <- list(locked = build(GENE_FRC, "locked"))
CFG <- CFG[!vapply(CFG, is.null, logical(1))]

## ---- 3. eight end points ------------------------------------------------------
res <- do.call(rbind, lapply(names(CFG), function(tag) {
  z <- CFG[[tag]]; E <- clin[match(z$id, clin$STUDY_NUMBER), ]
  do.call(rbind, lapply(names(ENDPT), function(ep) {
    sc <- ENDPT[[ep]][1]; tc <- ENDPT[[ep]][2]
    t <- as.numeric(E[[tc]]); s <- as.numeric(E[[sc]])
    ce <- Cidx(t, s, z$risk); h <- hr_sd(t, s, z$risk)
    data.frame(config = tag, n = length(z$id), n_genes = z$ngene, endpoint = ep,
               events = sum(s %in% 1), C = round(ce[1], 3),
               C_lo = round(ce[1] - 1.96*ce[2], 3), C_hi = round(ce[1] + 1.96*ce[2], 3),
               AUC3 = round(auc3(t, s, z$risk), 3),
               HR_perSD = round(h[1], 3), HR_lo = round(h[2], 3), HR_hi = round(h[3], 3),
               p = signif(h[4], 3))
  }))
}))

cat("\n===== locked five-gene score, replication cohort =====\n")
for (tag in names(CFG)) {
  cat("\n--- config:", tag, "(n =", CFG[[tag]]$n <- length(CFG[[tag]]$id),
      ", genes =", CFG[[tag]]$ngene, ") ---\n")
  print(res[res$config == tag, setdiff(names(res), c("config","n","n_genes"))], row.names = FALSE)
}
write.csv(res, file.path(resd, "mace_subtype_validation_3batch.csv"), row.names = FALSE)
main <- CFG[[1]]
scores <- data.frame(STUDY_NUMBER = main$id, batch = as.character(main$bat),
                     risk = main$risk, main$X[, SIG, drop = FALSE], check.names = FALSE)
write.csv(scores, file.path(resd, "external3batch_riskscore_and_genes.csv"), row.names = FALSE)
cat("\nwritten:\n -", file.path(resd, "mace_subtype_validation_3batch.csv"),
    "\n -", file.path(resd, "external3batch_riskscore_and_genes.csv"), "\n")
