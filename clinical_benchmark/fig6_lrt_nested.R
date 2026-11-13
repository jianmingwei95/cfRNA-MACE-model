## ============================================================================
## Figure 6D add-on | likelihood-ratio test for the increment of cfRNA over each clinical score
##   Nested Cox models in the replication cohort (n = 277): Surv ~ clinical score  vs
##   Surv ~ clinical score + locked cfRNA score (1 df). Clinical score = published-coefficient LP.
##   End points: ischaemic stroke (Fig 6D) and MACE (for the text).
## Output: 03_results/fig6/fig6_lrt_nested.csv
## ============================================================================
root <- "path/to/project/external"
src <- readLines(file.path(root, "clinical_benchmark", "benchmark_locked_weights.R"))
eval(parse(text = src[seq_len(grep("^## ---- 1\\) weights", src) - 1)]))   # helpers, rep, SR (clinical scores)
LR <- rbindlist(lapply(c("EP_ischemic_stroke", "EP_major"), function(ep) {
  o <- outcome(rep$STUDY_NUMBER, ep)
  rbindlist(lapply(names(SR), function(k) {
    d <- data.frame(t = o$t, e = o$e, clin = SR[[k]], cf = rep$risk)
    f0 <- coxph(Surv(t, e) ~ clin, data = d); f1 <- coxph(Surv(t, e) ~ clin + cf, data = d)
    a <- anova(f0, f1)
    data.table(endpoint = ep, score = k, n = nrow(d), events = sum(d$e), chisq = a$Chisq[2], df = a$Df[2], p_lrt = a$`Pr(>|Chi|)`[2],
               hr_cf = exp(coef(f1)[["cf"]]))
  })) }))
fwrite(LR, file.path(root, "03_results", "fig6", "fig6_lrt_nested.csv"))
print(LR[, .(endpoint, score, events, chisq = round(chisq, 2), p_lrt = signif(p_lrt, 3))])
