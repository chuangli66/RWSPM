
t0 = proc.time()[3]

rm(list=ls())
require(randtoolbox, quietly = T)
require(pracma, quietly = T)
require(clue, quietly=T)
require(energy, quietly = T)

library(randtoolbox)
library(pracma)
library(clue)
library(energy)
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
  return(randcovSTAT$statistic)
}

load("/.../ADNI_1_3.Rdata")
# load("/.../left_res_ori.Rdata")
load("/.../right_res_ori.Rdata")
###########
# Rdcov #
##########
### Load the statistic under the null distribution
quanenergydcov = read.table("./.../QuanenergyRankDcov-99999-1-15000.txt")
X <- ADNI_COMBINE[,-1]
rm(ADNI_COMBINE)
M = dim(quanenergydcov)[1]
o <- matrix(NA, nrow = 494465, ncol = 3)
for (l in seq(494465)) {
  ### global test
  # res = computestatisticrdcov(as.matrix(X[,l]),left_res_ori)
  res = computestatisticrdcov(as.matrix(X[,l]),right_res_ori)
  pval = length(which(quanenergydcov>res))/M
  o[l,] = c(l, res, pval)
}
# write.csv(o, paste0("/.../rdcov_global_left",l,".csv"), row.names = F)
write.csv(o, paste0("/.../rdcov_global_right",l,".csv"), row.names = F)



