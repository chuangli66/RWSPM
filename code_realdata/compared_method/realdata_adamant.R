###########
# AdaMant #
###########
library(adamant)

# load("./data/left_res_ori.Rdata")
load("./data/right_res_ori.Rdata")
load("./data/ADNI_1_3.Rdata")
X <- ADNI_COMBINE[,-1]
rm(ADNI_COMBINE)

### global test
adamant_pval = matrix(NA, nrow = 494465, ncol = 2)
for (l in seq(494465)) {
  colnames(adamant_pval) = c("snp.l","pval")
  # res = adamant(left_res_ori, X[,l], verbose = F,n_perms = 9999)
  res = adamant(right_res_ori, X[,l], verbose = F,n_perms = 9999)
  adamant_pval[l,] = c(l, res$P_val)
}
# write.csv(adamant_pval, paste0("/.../adamant_global_left.csv"), row.names = F)
write.csv(adamant_pval, paste0("/.../adamant_global_right.csv"), row.names = F)


