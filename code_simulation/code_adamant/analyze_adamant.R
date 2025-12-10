

t0 = proc.time()[3]

library(adamant)
## File loading path
root = "/.../Setting" 
## file storage path
save_root = paste0("/.../Setting",k1,"/")
## k1 is the number of Setting 1-7
k1=1
cat(paste0("set",k1),"\n")
dor = c()
for (l in 1:100) {
  cat(paste0("set",k1,"region:",l),"\n")
  load(paste0(root, k1, "/sim_", l , "/Y_slide",l,".RData"))
  load(paste0(root, k1, "/sim_", l , "/x.RData"))

  ###########
  # AdaMant #
  ###########
  n = nrow(x) # number of subjects (sample size): 200
  m = length(Y_slide)
  o.adamant = matrix(NA, nrow = m, ncol = 2)
  colnames(o.adamant) = c("region.j","pval")

  t1 = proc.time()[3]
  for(j in 1:m){
    if(j%%100 ==0){print(j)}
    res = adamant(Y_slide[[j]], x, n_perms = 199, verbose = F)
    if (is.na(res$P_val)){o.adamant[j,] = c(j,1)
    res$P_val = 1
    }else{o.adamant[j,] = c(j, res$P_val)}
    if(res$P_val == 0){res = adamant(Y_slide[[j]], x, n_perms = 10000, verbose = F)}
  }
  round(proc.time()[3] - t1)
  write.csv(o.adamant, paste0(save_root,"test_adamant",l,".csv"), row.names = F)
}

## global test
## k1 is the number of Setting 1-7
k1=1
cat(paste0("set",k1),"\n")
o.adamant = matrix(NA, nrow = 100, ncol = 2)
colnames(o.adamant) = c("times.j","pval")
for (l in seq(100)) {
  cat(paste0("set",k1,"times:",l),"\n")
  load(paste0(root, k1, "/sim_", l , "/y.RData"))
  load(paste0(root, k1, "/sim_", l , "/x.RData"))
  res = adamant(y, x, n_perms = 199, verbose = F)
  o.adamant[l,] = c(l, res$P_val)
  cat(paste0("set",k1),":",o.adamant[l,],"\n")
}
print(mean(o.adamant[,2] < 0.05))
write.csv(o.adamant, paste0(save_root,"adamant_global_set",k1,".csv"), row.names = F)


