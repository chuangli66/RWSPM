# require(randtoolbox, quietly = T)
# require(pracma, quietly = T)
# require(clue, quietly=T)
# require(energy, quietly = T)
library(pracma)
library(randtoolbox)
library(clue)
library(energy)
quanenergydcov=scan("/.../QuanenergyRankDcov-15000-1-9.txt")
M = length(quanenergydcov)


######### Computing Rank Energy Statistic #########
computestatisticrdcov=function(x,y,dim1=ncol(x),dim2=ncol(y),n=nrow(x),gridch=halton(n,dim1+dim2),gridch1=as.matrix(gridch[,(1:dim1)]),gridch2=as.matrix(gridch[,((dim1+1):(dim1+dim2))]))
{
  distmat1=matrix(0,nrow=n,ncol=n)
  for(i in 1:n)
    distmat1[i,]=apply((x[i,]-t(gridch1)),2,Norm)^2
  assignmentFUN1=solve_LSAP(distmat1)
  assignmentSOL1=cbind(seq_along(assignmentFUN1),assignmentFUN1)
  distmat2=matrix(0,nrow=n,ncol=n)
  for(i in 1:n)
    distmat2[i,]=apply((y[i,]-t(gridch2)),2,Norm)^2
  assignmentFUN2=solve_LSAP(distmat2)
  assignmentSOL2=cbind(seq_along(assignmentFUN2),assignmentFUN2)
  randcovSTAT=dcov.test(gridch1[assignmentSOL1[,2],],gridch2[assignmentSOL2[,2],], R=1)
  # randcovSTAT = dcov(gridch1[assignmentSOL1[,2],],gridch2[assignmentSOL2[,2],])
  return(randcovSTAT$statistic)
  # return(randcovSTAT)
}

## File loading path
root = "/.../Setting" 
## file storage path
save_root = paste0("/.../Setting",k1,"/")
## k1 is the number of Setting 1-7
k1=1
  cat(paste0("set",k1),"\n")
  t0 = proc.time()[3]
  for (l in 1:100) {
    cat(l,"\n")
    load(paste0(root, k1, "/sim_", l , "/Y_slide",l,".RData"))
    load(paste0(root, k1, "/sim_", l , "/x.RData"))

    n = nrow(x)
    ###########
    # rdcov #
    ###########

    m = length(Y_slide)
    o = matrix(NA, nrow = m, ncol = 3)
    colnames(o) = c("region.j","statistic","pval")

    # t1 = proc.time()[3]
    for(j in 1:m){
      t0 = Sys.time()
      res = computestatisticrdcov(x,Y_slide[[j]])
      pval = length(which(quanenergydcov>res))/M
      o[j,] = c(j,res, pval)
      Sys.time() - t0
      if(j%%100 == 0){cat(paste0(l,"-",j),"\n")}
    }
    write.csv(o, paste0(save_root,"test_rdcov",l,".csv"), row.names = F)
  }
  cat(paste0(round(proc.time()[3] - t0),"\n"))
# }


## global test
M = length(quanenergydcov)
## k1 is the number of Setting 1-7
k1=2
cat(paste0("set",k1),"\n")
o.rdcov = matrix(NA, nrow = 99, ncol = 3)
colnames(o.rdcov) = c("times.j","statistic","pval")
for (l in seq(100)) {
  cat(paste0("set",k1,"times:",l),"\n")
  load(paste0(root, k1, "/sim_", l , "/y.RData"))
  load(paste0(root, k1, "/sim_", l , "/x.RData"))
  t0 = Sys.time()
  res = computestatisticrdcov(x,as.matrix(y))
  pval = length(which(quanenergydcov>res))/M
  o.rdcov[l,] = c(l,res, pval)
  Sys.time() - t0
  cat(paste0("set",k1),":",o.rdcov[l,],"\n")
}
print(mean(o.rdcov[,3] < 0.05))
write.csv(o.rdcov, paste0(save_root,"rdcov_global_set",k1,".csv"), row.names = F)






