######### install packages for functions
# halton 
# install.packages("randtoolbox")
require(randtoolbox, quietly = T)

# Norm
# install.packages("pracma")
require(pracma, quietly = T)

# solve_LSAP
# install.packages("clue")
require(clue, quietly=T)

# dcov.test
# install.packages("energy")
require(energy, quietly = T)

######### Generating universal distribution ##########
gensamdistrdcov=function(N,dim1,dim2,niter=15000,fixgrid=halton(N,dim1+dim2),fixgrid1=as.matrix(fixgrid[,(1:dim1)]),fixgrid2=as.matrix(fixgrid[,((dim1+1):(dim1+dim2))]))
{
  tstat=numeric(niter)
  for(i in 1:niter)
  {
    ranper=sample(N)
    tstat[i]=dcov.test(fixgrid1,fixgrid2[ranper,],R=1)$statistic
    print(i)
  }
  return(tstat)
}
# load("E:/region_analysis/simulation/subr_idx_width27.RData")
dim1=1
dim2=9
quanenergydcov = gensamdistrdcov(200, dim1, dim2)
write.table(quanenergydcov, "/.../QuanenergyRankDcov-15000-1-9.txt",
            quote = F, col.names = F, row.names = F)

