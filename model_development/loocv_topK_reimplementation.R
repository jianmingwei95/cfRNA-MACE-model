#!/usr/bin/env Rscript
## Re-implementation of the discovery nested LOOCV (manuscript Methods) with a configurable number of genes K.
## Per held-out patient (344 iterations), using the remaining 343 only:
##   1) rank all genes by univariable-Cox C-index (direction-aware: max(C, 1 - C))
##   2) select the top K sequentially, skipping genes with |Pearson r| >= 0.7 to an already selected gene
##   3) z-score genes within the training fold, ridge Cox (glmnet alpha = 0; lambda by 10-fold CV within the training set)
##   4) linear predictor for the held-out patient
## Global (pooled) Harrell's C-index over the 344 out-of-fold scores, bootstrap 95% CI.
## Usage: Rscript loocv_topK_reimplementation.R <K>      (K = 10 reproduces the original; K = 5 is the new run)
suppressMessages({library(data.table);library(survival);library(glmnet);library(parallel)})
args<-commandArgs(TRUE); K<-if(length(args)) as.integer(args[1]) else 5L
CENTER<-if(length(args)>1) as.integer(args[2])==1 else TRUE   # FALSE = raw (uncentred) linear predictor, as in the original pipeline output
ROOT<-"path/to/project/external"
OUT<-file.path(ROOT,"03_results","loocv_topK"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
d<-fread("path/to/project/data/survival_data_cleanqc.csv")
time<-d$EP_major_time; ev<-d$EP_major
X<-as.matrix(d[,!c("STUDY_NUMBER","EP_major","EP_major_time"),with=FALSE]); n<-nrow(X)
cat(sprintf("n = %d, events = %d, genes = %d, K = %d\n",n,sum(ev),ncol(X),K))
## ---- vectorised leave-one-out univariable C-index for all genes and all folds ----
## Harrell's C over comparable pairs (i event, t_i < t_j); C_gene = (concordant + 0.5 ties) / pairs.
## Removing patient k only removes pairs that involve k, so C_{-k} = (N - n_k) / (D - d_k), exact.
ii<-which(ev==1); pr<-do.call(rbind,lapply(ii,function(i){j<-which(time>time[i] | (time==time[i] & ev==0));if(length(j)) cbind(i,j)}))
D<-nrow(pr); dk<-tabulate(c(pr[,1],pr[,2]),nbins=n)                      # pairs per patient
Nk<-matrix(0,n,ncol(X)); Ntot<-numeric(ncol(X))
chunks<-split(seq_len(ncol(X)),ceiling(seq_len(ncol(X))/500))
for(ch in chunks){A<-X[pr[,1],ch,drop=FALSE]; B<-X[pr[,2],ch,drop=FALSE]
  V<-(A>B)+0.5*(A==B)                                                    # pair x gene concordance
  Ntot[ch]<-colSums(V)
  Nk[,ch]<-rowsum(V,pr[,1],reorder=TRUE)[match(seq_len(n),sort(unique(pr[,1]))),,drop=FALSE] |> (\(m){m[is.na(m)]<-0;m})() +
           rowsum(V,pr[,2],reorder=TRUE)[match(seq_len(n),sort(unique(pr[,2]))),,drop=FALSE] |> (\(m){m[is.na(m)]<-0;m})()}
CLOO<-sweep(-Nk,2,Ntot,"+")/(D-dk)                                         # n x genes: C-index without patient k
SLOO<-pmax(CLOO,1-CLOO)
cat(sprintf("Vectorised LOO C-index matrix: %d x %d (comparable pairs %d)\n",nrow(SLOO),ncol(SLOO),D))
## correctness check against survival::concordance on a random fold and genes
set.seed(99); kk<-sample(n,1); gg<-sample(ncol(X),5); tr0<-setdiff(seq_len(n),kk)
chk<-sapply(gg,function(g)concordance(Surv(time[tr0],ev[tr0])~X[tr0,g],reverse=TRUE)$concordance)
cat(sprintf("Check fold %d: max |vectorised - survival::concordance| = %.2e\n",kk,max(abs(chk-CLOO[kk,gg]))))
one_fold<-function(k){tr<-setdiff(seq_len(n),k)
  ord<-order(-SLOO[k,])
  sel<-integer(0)
  for(g in ord){ if(length(sel)>=K) break
    if(length(sel)==0 || all(abs(cor(X[tr,g],X[tr,sel,drop=FALSE]))<0.7)) sel<-c(sel,g)}
  ## standardise within the training fold; apply training mean/SD to the held-out patient (as StandardScaler)
  mu<-colMeans(X[tr,sel,drop=FALSE]); sdv<-apply(X[tr,sel,drop=FALSE],2,sd)
  Ztr<-sweep(sweep(X[tr,sel,drop=FALSE],2,mu,"-"),2,sdv,"/"); Zk<-(X[k,sel]-mu)/sdv
  y<-Surv(time[tr],ev[tr]); set.seed(k)
  if(CENTER){cv<-cv.glmnet(Ztr,y,family="cox",alpha=0,nfolds=10,standardize=FALSE)
    lp<-as.numeric(predict(cv,matrix(Zk,nrow=1),s="lambda.min",type="link"))
  } else {cv<-cv.glmnet(X[tr,sel,drop=FALSE],y,family="cox",alpha=0,nfolds=10)
    lp<-as.numeric(predict(cv,X[k,sel,drop=FALSE],s="lambda.min",type="link"))}
  list(lp=lp,genes=colnames(X)[sel])}
t0<-Sys.time()
R<-mclapply(seq_len(n),one_fold,mc.cores=max(1,detectCores()-1))
cat(sprintf("LOOCV done in %.1f min\n",as.numeric(difftime(Sys.time(),t0,units="mins"))))
lp<-sapply(R,`[[`,"lp"); G<-sapply(R,function(r)paste(r$genes,collapse=","))
C<-concordance(Surv(time,ev)~lp,reverse=TRUE)$concordance
set.seed(1); bs<-replicate(1000,{i<-sample(n,replace=TRUE);concordance(Surv(time[i],ev[i])~lp[i],reverse=TRUE)$concordance})
cat(sprintf("[%s] Global LOOCV C-index (K = %d): %.3f (95%% CI %.3f-%.3f)\n",ifelse(CENTER,"centred","raw"),K,C,quantile(bs,.025),quantile(bs,.975)))
fq<-sort(table(unlist(lapply(R,`[[`,"genes"))),decreasing=TRUE)
cat("Selection frequency (top 15):\n");print(head(fq,15))
fwrite(data.table(K=K,C=C,lo=quantile(bs,.025),hi=quantile(bs,.975)),file.path(OUT,sprintf("loocv_top%d_%s_cindex.csv",K,ifelse(CENTER,"centred","raw"))))
fwrite(data.table(gene=names(fq),n_folds=as.integer(fq),pct=100*as.integer(fq)/n),file.path(OUT,sprintf("loocv_top%d_%s_selection_freq.csv",K,ifelse(CENTER,"centred","raw"))))
saveRDS(list(lp=lp,genes=G,K=K),file.path(OUT,sprintf("loocv_top%d_%s_oof.rds",K,ifelse(CENTER,"centred","raw"))))
