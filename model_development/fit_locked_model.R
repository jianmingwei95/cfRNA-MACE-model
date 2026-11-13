## ============================================================================
## Final five-gene model (SPINT2, EIF4E, TRIM25, EI24, F11R), discovery cohort (n = 344)
##   Ridge-penalised Cox model for MACE (glmnet, alpha = 0); penalty parameter chosen by
##   10-fold cross-validation (lambda.min). The fitted object, the discovery reference
##   distribution of the five genes (used for quantile mapping of the replication cohort)
##   and the tertile cut-offs of the discovery score are saved and used unchanged afterwards.
## Output: 03_results/disc344_to_external/model_5gene.rds
## ============================================================================
suppressMessages({ library(survival); library(glmnet) })
root <- "path/to/project/external"
resd <- file.path(root, "03_results", "disc344_to_external"); dir.create(resd, showWarnings = FALSE, recursive = TRUE)
SIG  <- c("SPINT2", "EIF4E", "TRIM25", "EI24", "F11R")
disc <- read.csv("path/to/project/data/survival_data_cleanqc.csv", check.names = FALSE)
Xtr  <- as.matrix(disc[, SIG]); ytr <- Surv(disc$EP_major_time, disc$EP_major)
set.seed(42)
cv   <- cv.glmnet(Xtr, ytr, family = "cox", alpha = 0, nfolds = 10)
fit  <- glmnet(Xtr, ytr, family = "cox", alpha = 0, lambda = cv$lambda.min)
lp   <- as.numeric(predict(fit, Xtr, type = "link"))
cat("lambda.min =", signif(cv$lambda.min, 5), "\n"); print(round(as.matrix(coef(fit)), 4))
saveRDS(list(sig = SIG, fit = fit, lambda = cv$lambda.min, cuts = quantile(lp, c(1/3, 2/3)), disc_ref = disc[, SIG, drop = FALSE]),
        file.path(resd, "model_5gene.rds"))
