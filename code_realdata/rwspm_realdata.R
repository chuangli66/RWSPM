## instal package RWSPM (/.../RWSPM_0.1.0.tar.gz)
library(RWSPM)
load("/.../ADNI_1_3.Rdata")
load("/.../left_res_G.Rdata")
# load("/.../right_res_G.Rdata")
X <- ADNI_COMBINE[,-1]
y_hipp = as.matrix(right_res_hipp)


cct.rwspm.pvalue.right <- matrix(NA,494465,2)
for (l in seq(494465)) {
  ## -----------------------------
  ## RegionWise statistical parametric mapping
  ## -----------------------------
  rw.pvalue = RWSPM(x = X[,l], y = y_hipp, M1 = 150, M2 = 100, window.width = 11)
  # write.csv(rw.pvalue$rw.pvalue, file = paste0("./.../left_rwspm_pvalue",l,".csv"),row.names = FALSE)
  write.csv(rw.pvalue$rw.pvalue, file = paste0("./.../right_rwspm_pvalue",l,".csv"),row.names = FALSE)
  ## -----------------------------
  ## Modified Truncated Cauchy Combination Test on the most significant cluster
  ## -----------------------------
  p.value <- RW.MTCCT(rw.pvalue$rw.pvalue, rw.pvalue$num.region.r, rw.pvalue$num.region.c, threshold = 10^-5)
  cct.rwspm.pvalue.right[i,] <- c(i,p.value)
}

# write.csv(cct.rwspm.pvalue.left, paste0("/.../adamant_global_left.csv"), row.names = F)
write.csv(cct.rwspm.pvalue.right, paste0("/.../adamant_global_right.csv"), row.names = F)


