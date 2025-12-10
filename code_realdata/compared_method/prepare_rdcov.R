# load("D:/lichuang/RWSPM/realdata/realdata_result_1215/left/rdcov/rdcov_pval_top1000.RData")
# ### 1 snps Y_1281*400 1.42 min
######### install packages for functions
# halton
install.packages("randtoolbox")
require(randtoolbox, quietly = T)

# Norm
install.packages("pracma")
require(pracma, quietly = T)

# solve_LSAP
install.packages("clue")
require(clue, quietly=T)

# dcov.test
install.packages("energy")
require(energy, quietly = T)

library(randtoolbox)
library(pracma)
library(clue)
library(energy)

######### Generating universal distribution ##########
library(snowfall)
gensamdistrdcov1 <- function(l, N, dim1, dim2) {
  fixgrid=halton(N, dim1+dim2)
  fixgrid1=as.matrix(fixgrid[, (1:dim1)])
  fixgrid2=as.matrix(fixgrid[, ((dim1+1):(dim1+dim2))])
    ranper = sample(N)
    tstat = dcov.test(fixgrid1, fixgrid2[ranper,], R=1)$statistic
  return(tstat)
}

sfInit(parallel = TRUE, cpus = 4)
sfLibrary(randtoolbox)
sfLibrary(energy)

dim1 = 1
dim2 = 15000
niter = 9999
N = 1270

t0 = Sys.time()
result <- sfLapply(1:niter, gensamdistrdcov1, N = N, dim1 = dim1, dim2 = dim2)
Sys.time() - t0
print(unlist(result))
quanenergydcov = unlist(result)
write.table(quanenergydcov, "./code_realdata/compared_method/QuanenergyRankDcov-99999-1-15000.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

sfStop()





