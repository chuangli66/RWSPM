library(GBJ)
## File loading path
root = "/.../Setting" 
## file storage path
save_root = paste0("/.../Setting",k1,"/")
## k1 is the number of Setting 1-7
k1=1
cat(paste0("set",k1),"\n")

for (l in 1:100) {
  cat(paste0("set",k1,"region:",l),"\n")
  load(paste0(root, k1, "/sim_", l , "/Y_slide",l,".RData"))
  load(paste0(root, k1, "/sim_", l , "/x.RData"))
  # Y = as.matrix(y); rm(y)
  n = nrow(x) # number of subjects (sample size): 200
  m = length(Y_slide)


  #######
  # GBJ #
  #######
  d = dim(Y_slide[[1]])[2]

  o.gbj = matrix(NA, nrow = m, ncol = (d+4))
  colnames(o.gbj) = c("region.j", paste0("glm.pval", 1:d),
                      "GBJ.stat", "GBJ.pvalue", "GBJ.message")

  in.GBJ = matrix(NA, ncol = 3, nrow = d)
  colnames(in.GBJ) = c("stat", "sigma", "pvalue")
  row.names(in.GBJ) = 1:d

  x = as.matrix(x)
  for(j in 1:m){
    o.gbj[j,1] = j
    Y = Y_slide[[j]]
    for(k in 1:d){
      y = Y[,k]
      res = glm(y ~ x)
      summ = summary(res)
      in.GBJ[k,] = c(summ$coefficients[2,3], # t-value
                     sqrt(var(res$residuals)), # sigma_j
                     summ$coefficients[2,4]) # p-value
    }
    o.gbj[j,2:(d+1)] = in.GBJ[,3]

    Y_adj = Y
    for(k in 1:d){
      Y_adj[,k] = Y[,k]/ in.GBJ[k,2]
    }
    cor.mat = cor(Y_adj)

    res = GBJ(in.GBJ[,1], cor_mat = cor.mat)

    o.gbj[j,(ncol(o.gbj)-2):ncol(o.gbj)] = c(res$GBJ, res$GBJ_pvalue, res$err_code)

  }
  print(c(l,cover_function_compare1(t(as.matrix(o.gbj[,"GBJ.pvalue"])))))
  oo.gbj = as.data.frame(o.gbj)s
  oo.gbj$num_signif = apply(o.gbj[,2:(d+1)],1,function(v) sum(as.numeric(as.character(v)) < 0.05))
  write.csv(oo.gbj, paste0(save_root,"test_GBJ",l,".csv"), row.names = F)
}
