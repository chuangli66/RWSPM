library(xtable)
library(igraph)

load("/.../snp_inf1_3.Rdata")

#############################
######## rwspm ############
#############################
pval_rwspm = read.csv("/.../rwspm_global_right.csv",header = TRUE)
top.rwspm.right <- cbind(snp_inf1,pval_rwspm[,2])
sort_index_rwspm.right = sort(pval_rwspm_top[,2],index.return = TRUE)
top_rwspm <- cbind(top.rwspm.right[sort_index_rwspm.right$ix[1:10],c(1,2,4,7)])
colnames(top_skpcr) = c("chr","snp","position","p.value")
write.csv(top_rwspm,file = "/.../right_hipp_rwspm.csv")

#############################
######## gbj ############
#############################
pval_gbj = read.csv("/.../gbj_global_right.csv",header = TRUE)
top.gbj.right <- cbind(snp_inf1,pval_gbj[,2])
sort_index_gbj.right = sort(pval_gbj_top[,2],index.return = TRUE)
top_gbj <- cbind(top.gbj.right[sort_index_gbj.right$ix[1:10],c(1,2,4,7)])
colnames(top_skpcr) = c("chr","snp","position","p.value")
write.csv(top_gbj,file = "/.../right_hipp_gbj.csv")

#############################
######## adamant ############
#############################
pval_adamant = read.csv("/.../adamant_global_right.csv",header = TRUE)
top.adamant.right <- cbind(snp_inf1,pval_adamant[,2])
sort_index_adamant.right = sort(pval_adamant_top[,2],index.return = TRUE)
top_adamant <- cbind(top.adamant.right[sort_index_adamant.right$ix[1:10],c(1,2,4,7)])
colnames(top_skpcr) = c("chr","snp","position","p.value")
write.csv(top_adamant,file = "/.../right_hipp_adamant.csv")

#############################
######## skpcr ##############
#############################
pval_skpcr = read.table("/.../skpcr_global_right.txt")
top.skpcr.right <- cbind(snp_inf1,pval_skpcr[,2])
sort_index_skpcr.right = sort(pval_skpcr[,2],index.return = TRUE)
top_skpcr <- cbind(top.skpcr.right[sort_index_skpcr.right$ix[1:10],c(1,2,4,7)])
colnames(top_skpcr) = c("chr","snp","position","p.value")
write.csv(top_skpcr,file = "/.../left_hipp_skpcr.csv")
#############################
######## fvgwas #############
#############################
pval_fvgwas = read.table("/.../fvgwas_global_right.txt")
top.fvgwas.right <- cbind(snp_inf1,pval_fvgwas[,2])
sort_index_fvgwas.right = sort(pval_fvgwas[,2],index.return = TRUE)
top_fvgwas <- cbind(top.fvgwas.right[sort_index_fvgwas.right$ix[1:10],c(1,2,4,7)])
colnames(top_fvgwas) = c("chr","snp","position","p.value")
write.csv(top_fvgwas,file = "/.../right_hipp_fvgwas.csv")
#############################
######## rdcov  #############
#############################
pval_fvgwas = read.table("/.../rdcov_global_right.txt")
top.rdcov.right <- cbind(snp_inf1,pval_rdcov[,3])
sort_index_rdcov.right = sort(pval_rdcov[,3],index.return = TRUE)
top_rdcov <- cbind(top.rdcov.right[sort_index_rdcov.right$ix[1:10],c(1,2,4,7)])
colnames(top_rdcov) = c("chr","snp","position","p.value")
write.csv(top_rdcov,file = "/.../right_hipp_rdcov.csv")
