library(RWSPM)
## File loading path
root = "/.../Setting" 
## file storage path
save_root = paste0("/.../Setting",k1,"/")
rw.mtcct = re
cat(paste0("set",k),"\n")
for (l in 1:100) {
  print(l)
  
  load(paste0(root, k, "sim_", l , "/y.RData"))
  load(paste0(root, k, "sim_", l , "/x.RData"))
  
  ## -----------------------------
  ## RegionWise statistical parametric mapping
  ## -----------------------------
  rw.pvalue = RWSPM(x = x, y = y, M1 = 150, M2 = 100, window.width = 20)
  head(rw.pvalue$rw.pvalue)
  write.csv(rw.pvalue$rw.pvalue, paste0(save_root,k,"/test_rwspm",l,".csv"), row.names = F)
  ## -----------------------------
  ## Modified Truncated Cauchy Combination Test on the most significant cluster
  ## -----------------------------
  rw.mtcct[l] <- RW.MTCCT(rw.pvalue$rw.pvalue, rw.pvalue$num.region.r, rw.pvalue$num.region.c, threshold = 0.05)
  
}
write.csv(rw.mtcct[l], paste0(save_root,k,"/rwspm_cluster_set",l,".csv"), row.names = F)

