## instal package RWSPM (/.../RWSPM_0.1.0.tar.gz)
library(RWSPM)
CCT <- function(pvals,weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }
  
  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }
  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }
  
  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5 - pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }
  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}
k1=2
root = paste0("./.../Setting",k1,"/")
### bcov
bcov.pvalue <- read.table(file = paste0(root,"/pval_bcov_set",k1,".csv") ,sep = ",",header=TRUE)
bcov_global = rep(1,100)
for (i in seq(100)) {
  bcov_global[i] = RW.MTCCT(as.matrix(bcov.pvalue[i,-1]),19,29)
}
### gbj
gbj.pvalue <- read.table(file = paste0(root,"/pval_gbj_set",k1,".csv"),sep = ",",header=TRUE)
gbj_global = rep(1,100)
for (i in seq(100)) {
  gbj_global[i] = CCT(as.matrix(gbj.pvalue[i,-1]))
}
write.csv(bcov_global, file = paste0(root,"rwspm_global_set",k1,".csv"))
write.csv(gbj_global, file = paste0(root,"gbj_global_set",k1,".csv"))
### adamant
adamant_global = read.csv(file = paste0(root,"adamant_global_set",k1,".csv"))
### rdcov
rdcov_global = read.csv(file = paste0(root,"rdcov_global_set",k1,".csv"))
### skpcr
skpcr_global = read.table(file = paste0(root,"skpcr_global_set",k1,".txt"))
### fvgwas
fvgwas_global = read.csv(file = paste0(root,"fvgwas_global_set",k1,".csv"), header = FALSE)

print(mean(bcov_global<0.05))
print(mean(gbj_global<0.05))
print(mean(adamant_global[,2] < 0.05))
print(mean(rdcov_global[,3] < 0.05))
print(mean(skpcr_global[,2] < 0.05))
print(mean(fvgwas_global < 0.05))