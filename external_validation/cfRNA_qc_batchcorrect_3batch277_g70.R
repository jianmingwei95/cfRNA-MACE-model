## cfRNA 三批次 QC + 批次校正 —— n=277, 基因过滤放宽到 >=70% (主分析口径)
##
## 与 cfRNA_qc_batchcorrect_3batch.R 唯一的区别:
##   在 cbind 合并之后、任何 QC 之前, 先剔除与 discovery 队列重叠的样本,
##   口径与 mace_subtype_validation_3batch.R / leakage_audit_3batch.R 完全一致。
##   (旧版 278 例里有 1 例与 discovery 重叠, 已剔除)
##
## 输出统一加 _3batch277g70 / combined3_277g70_ 后缀, 不覆盖任何既有结果。
## 另外重出 PCA 批次图(沿用 pca_batch_3batch_replot.R 的椭圆画法)与终点事件表。

suppressMessages({
  library(org.Hs.eg.db); library(AnnotationDbi)
  library(limma); library(sva); library(haven); library(dplyr)
})

root <- "path/to/project/external"
setwd(root)
RAW <- file.path(root, "01_raw_data")
OUT <- file.path(root, "02_processed_data")
FIG <- file.path(root, "04_figures")
RES <- file.path(root, "03_results")
DISC <- "path/to/project/data/survival_data_cleanqc.csv"

f1 <- file.path(RAW, "cfRNAs_MOK12069_combined_raw_counts.txt")   # 批次 A
f2 <- file.path(RAW, "cfRNAs_MOK12755_combined_raw_counts.txt")   # 批次 B
f3 <- file.path(RAW, "cfRNAs_MOK12754_combined_raw_counts.txt")   # 批次 C

## ---- 0. 读入三批原始 counts ----
c1 <- read.table(f1, header = TRUE, row.names = 1, check.names = FALSE)
c2 <- read.table(f2, header = TRUE, row.names = 1, check.names = FALSE)
c3 <- read.table(f3, header = TRUE, row.names = 1, check.names = FALSE)
cat("批次A MOK12069 raw:", dim(c1), "| 批次B MOK12755 raw:", dim(c2),
    "| 批次C MOK12754 raw:", dim(c3), "\n")

## ---- 1. 共同基因 (Ensembl) 取交集后 cbind ----
common <- Reduce(intersect, list(rownames(c1), rownames(c2), rownames(c3)))
cat("三批共同基因:", length(common), "\n")
combined <- cbind(c1[common, , drop = FALSE],
                  c2[common, , drop = FALSE],
                  c3[common, , drop = FALSE])
batch <- factor(c(rep("MOK12069", ncol(c1)), rep("MOK12755", ncol(c2)), rep("MOK12754", ncol(c3))),
                levels = c("MOK12069", "MOK12755", "MOK12754"))
## 列名统一成纯数字 STUDY_NUMBER (原脚本在第7步做, 这里提前, 以便按 ID 去重)
ids <- gsub("\\D+", "", colnames(combined))
stopifnot(sum(duplicated(ids)) == 0, sum(ids == "") == 0)   # 三批样本 ID 不得重叠
colnames(combined) <- ids; names(batch) <- ids
cat("合并矩阵:", dim(combined), "| A:", ncol(c1), "B:", ncol(c2), "C:", ncol(c3), "\n")

## ---- 1b. **先**剔除与 discovery 重复的样本, 之后所有 QC 都在干净集合上做 ----
disc_id <- as.integer(read.csv(DISC, check.names = FALSE)$STUDY_NUMBER)
ov <- intersect(as.integer(colnames(combined)), disc_id)
cat("与 discovery(", length(disc_id), ") 重复的样本:", length(ov),
    if (length(ov)) paste0(" -> ", paste(ov, collapse = ", ")) else "", " [在任何QC之前剔除]\n")
keep_ext <- !(as.integer(colnames(combined)) %in% ov)
combined <- combined[, keep_ext, drop = FALSE]; batch <- droplevels(batch[keep_ext])
cat("剔除后进入QC:", ncol(combined), "|",
    paste(names(table(batch)), table(batch), collapse = " "), "\n")

## ---- 2. 删除 ERCC spike-in ----
ercc <- grep("^ERCC", rownames(combined), value = TRUE)
combined <- combined[setdiff(rownames(combined), ercc), ]
cat("去 ERCC 后:", dim(combined), "| ERCC:", length(ercc), "\n")

## ---- 3. 样本 QC: 低于总reads第5百分位 剔除, 且检出基因数 > 10000 ----
total_reads <- colSums(combined)
thr <- quantile(total_reads, 0.05)
low <- which(total_reads < thr)
kept <- if (length(low) > 0) setdiff(seq_len(ncol(combined)), low) else seq_len(ncol(combined))
combined <- combined[, kept]; batch <- batch[kept]
keep2 <- colSums(combined > 0) > 10000
combined <- combined[, keep2]; batch <- batch[keep2]
cat("样本QC后:", dim(combined), "| 5%阈值=", round(thr),
    "| 各批保留:", paste(names(table(batch)), table(batch), collapse = " "), "\n")

## ---- 4. 基因 QC: CPM > 0.5 且在 >= 75% 样本中满足 ----
total_counts <- colSums(combined)                 # 样本文库大小 (CPM 分母, 沿用原脚本)
cpm <- sweep(combined, 2, total_counts, "/") * 1e6
GENE_FRC <- 0.70                                  # 主分析口径 (75% 下 EI24 通不过)
gp  <- rowSums(cpm > 0.5) >= GENE_FRC * ncol(combined)
filtered_genes <- combined[gp, ]
cat(sprintf("基因CPM过滤后 (>=%.0f%%): %d %d\n", 100*GENE_FRC, nrow(filtered_genes), ncol(filtered_genes)))

## ---- 5. 注释: 仅 protein-coding, 去重 symbol ----
filtered_genes$symbol   <- mapIds(org.Hs.eg.db, keys=rownames(filtered_genes), column="SYMBOL",   keytype="ENSEMBL", multiVals="first")
filtered_genes$genetype <- mapIds(org.Hs.eg.db, keys=rownames(filtered_genes), column="GENETYPE", keytype="ENSEMBL", multiVals="first")
filtered_genes <- filtered_genes[filtered_genes$genetype %in% "protein-coding", ]
filtered_genes <- filtered_genes[complete.cases(filtered_genes$symbol), ]
filtered_genes <- filtered_genes[!duplicated(filtered_genes$symbol), ]
rownames(filtered_genes) <- filtered_genes$symbol
filtered_genes$symbol <- NULL; filtered_genes$genetype <- NULL
cat("protein-coding+去重后:", dim(filtered_genes), "\n")

## ---- 6. 去核糖体基因 (RPL*/RPS*) ----
ribo <- c(grep("^RPL", keys(org.Hs.eg.db, keytype="SYMBOL"), value=TRUE),
          grep("^RPS", keys(org.Hs.eg.db, keytype="SYMBOL"), value=TRUE))
qc <- filtered_genes[!rownames(filtered_genes) %in% ribo, ]
cat("去核糖体后 (最终QC矩阵 genes x samples):", dim(qc), "| ribo:", nrow(filtered_genes)-nrow(qc), "\n")

## ---- 7. 落盘 QC 矩阵 (列名已是 STUDY_NUMBER) ----
stopifnot(identical(colnames(qc), names(batch)))
saveRDS(qc, file.path(OUT, "filtered_genes_no_ribo_3batch277g70.rds"))
write.csv(qc, file.path(OUT, "filtered_genes_no_ribo_3batch277g70.csv"))

## ---- 8. 归一化: CPM -> 分位数 -> log2(x+1) ----
cpm_qc <- sweep(qc, 2, total_counts, "/") * 1e6
qn     <- as.data.frame(normalizeQuantiles(cpm_qc))
logQ   <- as.data.frame(t(apply(qn, 1, function(x) log2(x + 1))))
colnames(logQ) <- colnames(qn)
saveRDS(logQ, file.path(OUT, "combined3_277g70_logQ.rds"))
write.csv(logQ, file.path(OUT, "combined3_277g70_logQ.csv"))

## ---- 9. 批次校正 ----
mat <- as.matrix(logQ)
ok <- apply(mat, 1, function(g) all(tapply(g, batch, function(v) sd(v) > 0)))
cat("批次内零方差基因剔除:", sum(!ok), "| 进入校正基因:", sum(ok), "\n")
mat <- mat[ok, ]

# 主方案: ComBat (仅按批次, 不放入结局, 避免泄漏)
combat_mat <- ComBat(dat = mat, batch = batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
logQ_combat <- as.data.frame(combat_mat)
saveRDS(logQ_combat, file.path(OUT, "combined3_277g70_logQ_combat.rds"))
write.csv(logQ_combat, file.path(OUT, "combined3_277g70_logQ_combat.csv"))

# 备选方案: limma::removeBatchEffect
logQ_rbe <- as.data.frame(removeBatchEffect(mat, batch = batch))
saveRDS(logQ_rbe, file.path(OUT, "combined3_277g70_logQ_removeBatchEffect.rds"))
write.csv(logQ_rbe, file.path(OUT, "combined3_277g70_logQ_removeBatchEffect.csv"))

## ---- 10. 批次元数据表 ----
meta <- data.frame(STUDY_NUMBER = colnames(qc), batch = as.character(batch),
                   lib_size = total_counts, genes_detected = colSums(qc > 0))
write.csv(meta, file.path(OUT, "combined3_277g70_batch_metadata.csv"), row.names = FALSE)

## ---- 11. PCA 校正前/后验证 (椭圆版, 同 pca_batch_3batch_replot.R) ----
pca_pc <- function(m) {
  p <- prcomp(t(m), scale. = FALSE)
  list(x = p$x, ve = round(100 * (p$sdev^2)/sum(p$sdev^2), 1))
}
pcaB <- pca_pc(mat); pcaA <- pca_pc(combat_mat)
r2 <- function(pc) summary(lm(pc ~ batch))$r.squared

pal  <- c(MOK12069 = "#1F4E79", MOK12755 = "#C0504D", MOK12754 = "#4E8542")
cols <- pal[as.character(batch)]
ell <- function(ctr, S, rad, n = 200) {
  th <- seq(0, 2*pi, length.out = n)
  t(ctr + rad * (chol(S) |> t()) %*% rbind(cos(th), sin(th)))
}
rad95 <- sqrt(qchisq(0.95, df = 2))
cen_sep <- function(X) {
  ctrs <- t(sapply(levels(batch), function(b) colMeans(X[batch == b, 1:2, drop = FALSE])))
  mean(dist(ctrs))
}
panel <- function(p, ttl) {
  X <- p$x[, 1:2, drop = FALSE]
  polys <- lapply(levels(batch), function(b) {
    Xb <- X[batch == b, , drop = FALSE]; ell(colMeans(Xb), cov(Xb), rad95)
  })
  xl <- range(c(X[,1], sapply(polys, function(q) q[,1])))
  yl <- range(c(X[,2], sapply(polys, function(q) q[,2])))
  plot(X, type = "n", axes = FALSE, ann = FALSE, xlim = xl, ylim = yl)
  usr <- par("usr"); rect(usr[1], usr[3], usr[2], usr[4], col = "grey97", border = NA)
  grid(col = "white", lty = 1, lwd = 1.1)
  points(X, col = adjustcolor(cols, alpha.f = 0.55), pch = 19, cex = 0.8)
  for (b in levels(batch)) {
    Xb <- X[batch == b, , drop = FALSE]; ctr <- colMeans(Xb); S <- cov(Xb); nb <- nrow(Xb)
    polygon(ell(ctr, S, rad95), border = pal[b], lwd = 2, lty = 2)
    polygon(ell(ctr, S/nb, rad95), border = pal[b], lwd = 2.2,
            col = adjustcolor(pal[b], alpha.f = 0.30))
    points(ctr[1], ctr[2], pch = 4, lwd = 2.4, cex = 1.15, col = pal[b])
  }
  axis(1, lwd = 0, lwd.ticks = 1.1, cex.axis = 0.9)
  axis(2, lwd = 0, lwd.ticks = 1.1, las = 1, cex.axis = 0.9)
  box(bty = "o", lwd = 1.6)
  title(xlab = sprintf("PC1 (%.1f%%)", p$ve[1]),
        ylab = sprintf("PC2 (%.1f%%)", p$ve[2]), line = 2.3)
  title(main = ttl, line = 1.9, cex.main = 1.05)
  mtext(sprintf("PC1~batch R2=%.3f   PC2~batch R2=%.3f   centroid sep=%.1f",
                r2(p$x[,1]), r2(p$x[,2]), cen_sep(p$x)),
        side = 3, line = 0.5, cex = 0.78, col = "grey25")
  legend("topright", legend = levels(batch), col = pal[levels(batch)],
         pch = 19, bty = "o", box.lwd = 1, bg = "white", cex = 0.82, pt.cex = 0.95)
}
for (d in list(c(file.path(FIG, "pca_batch_3batch277g70.pdf"), "pdf"),
               c(file.path(FIG, "pca_batch_3batch277g70.png"), "png"))) {
  if (d[2] == "pdf") pdf(d[1], width = 10.5, height = 5.2) else png(d[1], width = 1900, height = 950, res = 170)
  par(mfrow = c(1,2), mar = c(4.2, 4.4, 4, 1.2), mgp = c(3, 0.6, 0), tcl = -0.3)
  panel(pcaB, "Before correction"); panel(pcaA, "After ComBat")
  dev.off()
}
cat(sprintf("PC1~batch R^2  校正前=%.3f -> ComBat后=%.3f\n", r2(pcaB$x[,1]), r2(pcaA$x[,1])))
cat(sprintf("PC2~batch R^2  校正前=%.3f -> ComBat后=%.3f\n", r2(pcaB$x[,2]), r2(pcaA$x[,2])))
cat(sprintf("质心平均两两距离  校正前=%.2f -> 校正后=%.2f\n", cen_sep(pcaB$x), cen_sep(pcaA$x)))

## ---- 12. 合并临床生存终点 + 终点事件表 ----
sav <- file.path(RAW, "data20240320.sav")
if (file.exists(sav)) {
  clinical <- read_sav(sav); clinical$STUDY_NUMBER <- as.integer(clinical$STUDY_NUMBER)
  nums <- as.integer(colnames(qc))
  cat("== ID 匹配: cfRNA", length(nums), "| 匹配临床", sum(nums %in% clinical$STUDY_NUMBER),
      "| 未匹配", sum(!nums %in% clinical$STUDY_NUMBER), "\n")
  stopifnot(!any(nums %in% disc_id))   # 再次确认与 discovery 不相交

  surv <- clinical %>% dplyr::select(STUDY_NUMBER, EP_major, EP_major_time) %>%
    filter(STUDY_NUMBER %in% nums)
  bt <- data.frame(STUDY_NUMBER = nums, batch = as.character(batch))
  m  <- merge(surv, bt, by = "STUDY_NUMBER")
  cat("== batch x EP_major 事件平衡 ==\n"); print(table(m$batch, m$EP_major, useNA = "ifany"))

  tp <- as.data.frame(t(logQ_combat)); tp$STUDY_NUMBER <- as.integer(colnames(logQ_combat))
  sd_ <- merge(tp, surv, by = "STUDY_NUMBER") %>% filter(!is.na(EP_major) & !is.na(EP_major_time))
  sd_ <- merge(sd_, bt, by = "STUDY_NUMBER")
  saveRDS(sd_, file.path(OUT, "surivival_dataqc_3batch277g70.rds"))
  write.csv(sd_, file.path(OUT, "surivival_dataqc_3batch277g70.csv"), row.names = FALSE)
  cat("生存表(ComBat校正) surivival_dataqc_3batch277g70 dim:", dim(sd_),
      "| 事件:", sum(sd_$EP_major == 1), "| 删失:", sum(sd_$EP_major == 0), "\n")

  ## 终点事件表 (对应旧的 endpoint_counts_3batch278*.csv)
  ENDPT <- list(
    `Ischemic stroke`       = c("EP_ischemic_stroke", "EP_ischemic_stroke_time"),
    `Stroke (any)`          = c("EP_stroke",          "EP_stroke_time"),
    `Myocardial infarction` = c("EP_MI",              "EP_MI_time"),
    `CV death`              = c("EP_CVdeath",         "EP_CVdeath_time"),
    `MACE (EP_major)`       = c("EP_major",           "EP_major_time"),
    `MACE (EP_composite)`   = c("EP_composite",       "EP_composite_time"),
    `CAD`                   = c("EP_CAD",             "EP_CAD_time"),
    `Peripheral`            = c("EP_peripheral",      "EP_peripheral_time"),
    `Coronary`              = c("EP_coronary",        "EP_coronary_time"))
  Eclin <- clinical[match(nums, clinical$STUDY_NUMBER), ]
  bvec  <- as.character(batch)
  tab <- do.call(rbind, lapply(names(ENDPT), function(ep) {
    s <- as.numeric(Eclin[[ENDPT[[ep]][1]]]); t <- as.numeric(Eclin[[ENDPT[[ep]][2]]])
    ev <- s %in% 1
    data.frame(Endpoint = ep, Var = ENDPT[[ep]][1], N = length(nums),
               Events = sum(ev), Censored = sum(s %in% 0),
               Rate = sprintf("%.1f%%", 100*sum(ev)/length(nums)),
               MedianFU_years = round(median(t, na.rm = TRUE), 2),
               MaxFU_years    = round(max(t, na.rm = TRUE), 2),
               MOK12069 = sum(ev & bvec == "MOK12069"),
               MOK12755 = sum(ev & bvec == "MOK12755"),
               MOK12754 = sum(ev & bvec == "MOK12754"),
               check.names = FALSE)
  }))
  write.csv(tab, file.path(RES, "endpoint_counts_3batch277g70.csv"), row.names = FALSE)
  write.csv(tab[, c("Endpoint","MOK12069","MOK12755","MOK12754")] |>
              transform(Total = tab$Events),
            file.path(RES, "endpoint_counts_3batch277g70_bybatch.csv"), row.names = FALSE)
  cat("== 终点事件表 (n =", length(nums), ") ==\n"); print(tab, row.names = FALSE)
}

cat("\n完成. 关键输出:\n",
    " 02_processed_data/filtered_genes_no_ribo_3batch277g70.csv/.rds\n",
    " 02_processed_data/combined3_277g70_logQ{,_combat,_removeBatchEffect}.csv/.rds\n",
    " 02_processed_data/combined3_277g70_batch_metadata.csv\n",
    " 02_processed_data/surivival_dataqc_3batch277g70.csv/.rds\n",
    " 03_results/endpoint_counts_3batch277g70{,_bybatch}.csv\n",
    " 04_figures/pca_batch_3batch277g70.pdf/png\n")
