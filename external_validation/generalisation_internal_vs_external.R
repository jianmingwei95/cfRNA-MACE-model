## ============================================================================
## 5 基因 signature 的泛化落差: 内部 (discovery 344) vs 外部 (277), 8 个终点
##   1 Apparent            discovery 全体训练、同一批人评估 (乐观上界)
##   2 Cross-validated     discovery 10 折样本外预测 (诚实的内部估计)
##   3 External A          外部 277, 冻结系数 + 转导式标准化 (当前主结果)
##   4 External C          外部 277, 冻结系数 + 单样本冻结标准化 (前瞻式下界)
## 全部为同一个 5 基因模型 (alpha=0.1, 训练结局 EP_major), 系数除交叉验证外从不重估
## ============================================================================
root <- "path/to/project/external"
OUT <- file.path(root,"03_results"); resd <- file.path(OUT,"disc344_to_external")
I <- read.csv(file.path(OUT,"internal_discovery_5gene_endpoints.csv"), check.names=FALSE)
X <- read.csv(file.path(resd,"frozen_normalization_sensitivity_70pct.csv"), check.names=FALSE)
MAP <- c(`Ischaemic stroke`="Ischemic stroke", `Stroke (any)`="Stroke (any)",
         `Coronary artery disease`="CAD", `Coronary events`="Coronary",
         `Myocardial infarction`="Myocardial infarction",
         `MACE (training endpoint)`="MACE (EP_major)", `Cardiovascular death`="CV death",
         `Peripheral events`="Peripheral")
g <- function(sch, ep, col) { r <- X[X$scheme==sch & X$endpoint==MAP[[ep]], ]; if (!nrow(r)) NA else r[[col]] }
RES <- do.call(rbind, lapply(I$endpoint, function(ep) {
  i <- I[I$endpoint==ep, ]
  data.frame(endpoint=ep,
    ev_int=i$events, C_app=i$C_app, HR_app=i$HR_app, p_app=i$p_app,
    C_cv=i$C_cv, HR_cv=i$HR_cv, p_cv=i$p_cv,
    ev_ext=g("A transductive (main)",ep,"events"),
    C_extA=g("A transductive (main)",ep,"C"), HR_extA=g("A transductive (main)",ep,"HR_perSD"),
    p_extA=g("A transductive (main)",ep,"p"),
    C_extC=g("C per-sample frozen (prospective)",ep,"C"), HR_extC=g("C per-sample frozen (prospective)",ep,"HR_perSD"),
    p_extC=g("C per-sample frozen (prospective)",ep,"p"))
}))
RES$optimism  <- round(RES$C_app - RES$C_cv, 3)      # 内部乐观偏倚
RES$gen_gap   <- round(RES$C_cv  - RES$C_extC, 3)    # 内部诚实估计 -> 外部前瞻估计
RES <- RES[order(-RES$C_cv), ]
print(RES[,c("endpoint","ev_int","C_app","C_cv","optimism","ev_ext","C_extA","C_extC","gen_gap")],
      row.names=FALSE)
write.csv(RES, file.path(OUT,"generalisation_internal_vs_external.csv"), row.names=FALSE)
cat("\n== 汇总 ==\n")
cat(sprintf("内部乐观偏倚 (C_app - C_cv):   中位 %.3f, 范围 %.3f ~ %.3f\n",
    median(RES$optimism), min(RES$optimism), max(RES$optimism)))
cat(sprintf("泛化落差 (C_cv - C_extC):      中位 %.3f, 范围 %.3f ~ %.3f\n",
    median(RES$gen_gap), min(RES$gen_gap), max(RES$gen_gap)))
cat("\nwritten:", file.path(OUT,"generalisation_internal_vs_external.csv"), "\n")
