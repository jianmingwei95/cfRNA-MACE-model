## ============================================================================
## 5 基因 signature 在内部训练集 (discovery n = 344) 上、8 个终点的表现
##   模型: Cox glmnet alpha=0.1, lambda = cv.glmnet 10 折 lambda.min, 训练结局 EP_major
##   两套指标:
##     Apparent          —— 用全体 344 训练、又在同一批人上评估 (乐观, 会高估)
##     Cross-validated   —— 10 折: 每折在 9/10 上重训 (同样以 EP_major 为结局),
##                          在留出折上打分, 汇总后再对各终点评估 (去掉乐观偏倚)
## 输出: 03_results/internal_discovery_5gene_endpoints.csv + xlsx
## ============================================================================
suppressMessages({ library(haven); library(survival); library(glmnet); library(survivalROC) })
set.seed(42)
root <- "path/to/project/external"
OUT <- file.path(root,"03_results")
SIG <- c("SPINT2","EIF4E","TRIM25","EI24","F11R"); ALPHA <- 0.1
num <- function(x) as.numeric(unclass(x))

disc <- read.csv("path/to/project/data/survival_data_cleanqc.csv", check.names=FALSE)
cl <- read_sav(file.path(root,"01_raw_data","data20240320.sav")); cl$STUDY_NUMBER <- as.integer(cl$STUDY_NUMBER)
D <- cl[match(as.integer(disc$STUDY_NUMBER), cl$STUDY_NUMBER), ]
N <- nrow(disc); cat("discovery n =", N, "\n")

EPS <- list(`Ischaemic stroke`=c("EP_ischemic_stroke","EP_ischemic_stroke_time"),
            `Stroke (any)`=c("EP_stroke","EP_stroke_time"),
            `Coronary artery disease`=c("EP_CAD","EP_CAD_time"),
            `Coronary events`=c("EP_coronary","EP_coronary_time"),
            `Myocardial infarction`=c("EP_MI","EP_MI_time"),
            `MACE (training endpoint)`=c("EP_major","EP_major_time"),
            `Cardiovascular death`=c("EP_CVdeath","EP_CVdeath_time"),
            `Peripheral events`=c("EP_peripheral","EP_peripheral_time"))

X <- as.matrix(disc[, SIG, drop=FALSE]); y <- Surv(disc$EP_major_time, disc$EP_major)
cv <- cv.glmnet(X, y, family="cox", alpha=ALPHA, nfolds=10)
fit <- glmnet(X, y, family="cox", alpha=ALPHA, lambda=cv$lambda.min)
cat("lambda.min =", signif(cv$lambda.min,5), "| 系数:\n")
print(round(as.matrix(predict(fit, type="coefficients", s=cv$lambda.min))[SIG,1], 4))
risk_app <- as.numeric(predict(fit, X, type="link"))

## ---- 10 折交叉验证的样本外风险分 (训练结局始终是 EP_major) ----
folds <- sample(rep(1:10, length.out=N))
risk_cv <- numeric(N)
for (f in 1:10) {
  tr <- folds != f; te <- !tr
  cvf <- cv.glmnet(X[tr,,drop=FALSE], y[tr], family="cox", alpha=ALPHA, nfolds=10)
  ff  <- glmnet(X[tr,,drop=FALSE], y[tr], family="cox", alpha=ALPHA, lambda=cvf$lambda.min)
  risk_cv[te] <- as.numeric(predict(ff, X[te,,drop=FALSE], type="link"))
}
cat("apparent vs CV 风险分 Spearman =", round(cor(risk_app, risk_cv, method="spearman"),3), "\n")

metr <- function(t, s, r) {
  ok <- !is.na(t) & !is.na(s) & !is.na(r)
  if (sum(s[ok] %in% 1) < 3) return(rep(NA, 8))
  m <- summary(coxph(Surv(t[ok], s[ok]) ~ r[ok]))
  ms <- summary(coxph(Surv(t[ok], s[ok]) ~ scale(r[ok])[,1]))
  a <- tryCatch(survivalROC(t[ok], s[ok], r[ok], predict.time=3, method="KM")$AUC, error=function(e) NA)
  c(m$concordance[1], m$concordance[1]-1.96*m$concordance[2], m$concordance[1]+1.96*m$concordance[2],
    a, ms$conf.int[1,1], ms$conf.int[1,3], ms$conf.int[1,4], ms$coefficients[1,5])
}
RES <- do.call(rbind, lapply(names(EPS), function(ep) {
  t <- num(D[[EPS[[ep]][2]]]); s <- num(D[[EPS[[ep]][1]]])
  a <- metr(t, s, risk_app); v <- metr(t, s, risk_cv)
  data.frame(endpoint=ep, n=N, events=sum(s %in% 1), rate=100*sum(s %in% 1)/N,
             C_app=round(a[1],3), C_app_lo=round(a[2],3), C_app_hi=round(a[3],3),
             AUC3_app=round(a[4],3), HR_app=round(a[5],3), HR_app_lo=round(a[6],3),
             HR_app_hi=round(a[7],3), p_app=signif(a[8],3),
             C_cv=round(v[1],3), C_cv_lo=round(v[2],3), C_cv_hi=round(v[3],3),
             AUC3_cv=round(v[4],3), HR_cv=round(v[5],3), HR_cv_lo=round(v[6],3),
             HR_cv_hi=round(v[7],3), p_cv=signif(v[8],3))
}))
RES <- RES[order(-RES$C_app), ]
print(RES[,c("endpoint","events","C_app","AUC3_app","HR_app","p_app","C_cv","HR_cv","p_cv")], row.names=FALSE)
write.csv(RES, file.path(OUT,"internal_discovery_5gene_endpoints.csv"), row.names=FALSE)
cat("\nwritten:", file.path(OUT,"internal_discovery_5gene_endpoints.csv"), "\n")
