## instal package RWSPM (/.../RWSPM_0.1.0.tar.gz)
library(RWSPM)
k1=2
save_root = "/.../Setting"
### bcov
bcov_pvalue_mat <- NULL
for (i in seq(100)) {
  bcov_pvalue = read.csv(paste0(save_root,k1,"/test_BCov",i,".csv"),head = TRUE)
  bcov_pvalue_mat = rbind(bcov_pvalue_mat,bcov_pvalue[,3])
}
write.csv(bcov_pvalue_mat, file = paste0(save_root,k1,"/pval_bcov_set",k1,".csv"))

### adamant
adamant_pvalue_mat <- NULL
for (i in seq(100)) {
  adamant_pvalue = read.csv(paste0(save_root,k1,"/test_adamant",i,".csv"),head = TRUE)
  adamant_pvalue_mat = rbind(adamant_pvalue_mat,adamant_pvalue[,2])
}
write.csv(adamant_pvalue_mat, file = paste0(save_root,k1,"/pval_adamant_set",k1,".csv"))

### GBJ
GBJ_pvalue_mat <- NULL
for (i in seq(100)) {
  file_path = paste0(save_root,k1,"/test_GBJ",i,".csv")
  if(file.exists(file_path)){
    GBJ_pvalue = read.csv(file_path,head = TRUE)
    GBJ_pvalue_mat = rbind(GBJ_pvalue_mat,GBJ_pvalue[,"GBJ.pvalue"])
  }else{print("file doesn't exist")}
  print(i)
}
write.csv(GBJ_pvalue_mat, file = paste0(save_root,k1,"/pval_gbj_set",k1,".csv"))

### skpcr
skpcr_pvalue_mat <- NULL
for (i in seq(100)) {
  skpcr_pvalue = read.table(paste0(save_root,k1,"/test_skpcr",i,".txt"))
  skpcr_pvalue_mat = rbind(skpcr_pvalue_mat,skpcr_pvalue[,2])
}
write.csv(skpcr_pvalue_mat, file = paste0(save_root,k1,"/pval_skpcr_set",k1,".csv"))

### rdcov
rdcov_pvalue_mat <- NULL
for (i in seq(100)) {
  rdcov_pvalue = read.csv(paste0(save_root,k1,"/test_rdcov",i,".csv"),head = TRUE)
  rdcov_pvalue_mat = rbind(rdcov_pvalue_mat,rdcov_pvalue[,3])
}
write.csv(rdcov_pvalue_mat, file = paste0(save_root,k1,"/pval_rdcov_set",k1,".csv"))

### fvgwas
fvgwas_pvalue_mat <- NULL
for (i in seq(100)) {
  fvgwas_pvalue = read.csv(paste0(save_root,k1,"/test_fvgwas",i,".csv"),header = FALSE)
  fvgwas_pvalue_mat = rbind(fvgwas_pvalue_mat,fvgwas_pvalue[,2])
}
write.csv(fvgwas_pvalue_mat, file = paste0(save_root,k1,"/pval_fvgwas_set",k1,".csv"))
