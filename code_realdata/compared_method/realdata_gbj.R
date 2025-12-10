rm(list=ls())
library(GBJ)
library(RWSPM)
load("/.../left_res_hipp_width3.RData")
load("/.../right_res_hipp_width3.RData")
load("/.../ADNI_1_3.Rdata")
X <- ADNI_COMBINE[,-1]
rm(ADNI_COMBINE)
CCT <- function(pvals,weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }
  
  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }
  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }
  
  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5 - pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }
  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}

#######
# GBJ #
#######
gbj_function <- function(l) {
    library(GBJ)
    n = nrow(X) # number of subjects (sample size): 200
    m = length(left_res_hipp_width3)

    #######
    # GBJ #
    #######
    d = dim(left_res_hipp_width3[[1]])[2]

    o.gbj = matrix(NA, nrow = m, ncol = (d + 4))
    colnames(o.gbj) = c("region.j", paste0("glm.pval", 1:d),
                        "GBJ.stat", "GBJ.pvalue", "GBJ.message")

    in.GBJ = matrix(NA, ncol = 3, nrow = d)
    colnames(in.GBJ) = c("stat", "sigma", "pvalue")
    row.names(in.GBJ) = 1:d

    x = as.matrix(X[,gbj.sort$ix[l]])
    for (j in 1:m) {
      if (j %% 100 == 0) { print(j) }
      o.gbj[j, 1] = j
      Y = left_res_hipp_width3[[j]]
      t0 = Sys.time()
      for (k in 1:d) {
        y = Y[, k]
        res = glm(y ~ x)
        summ = summary(res)
        in.GBJ[k, ] = c(summ$coefficients[2, 3], # t-value
                        sqrt(var(res$residuals)), # sigma_j
                        summ$coefficients[2, 4]) # p-value
      }
      print(Sys.time() - t0)
      #Time difference of 1.731338 secs
      o.gbj[j, 2:(d + 1)] = in.GBJ[, 3]

      Y_adj = Y
      for (k in 1:d) {
        Y_adj[, k] = Y[, k] / in.GBJ[k, 2]
      }

      cor.mat = cor(Y_adj)
      t0 = Sys.time()
      res = GBJ(in.GBJ[, 1], cor_mat = cor.mat)
      print(Sys.time() - t0)
      #Time difference of 4.466828 mins
      o.gbj[j, (ncol(o.gbj) - 2):ncol(o.gbj)] = c(res$GBJ, res$GBJ_pvalue, res$err_code)
    }

    oo.gbj = as.data.frame(o.gbj)
    oo.gbj$num_signif = apply(o.gbj[, 2:(d + 1)], 1, function(v) sum(as.numeric(as.character(v)) < 0.05))

    # write.csv(oo.gbj, paste0("/.../left_gbj_pvalue",l,".csv"), row.names = F)
    write.csv(oo.gbj, paste0("/.../right_gbj_pvalue",l,".csv"), row.names = F)
    return(oo.gbj)
}

for (l in seq(494465)) {
  gbj_function(l)
}

# gbj.pvalue.list.left <- list()
gbj.pvalue.list.right <- list()
# cct.gbj.pvalue.left <- matrix(NA,494465,2)
cct.gbj.pvalue.right <- matrix(NA,494465,2)
for (i in seq(1,494465)) {
  # gbj.pvalue.load.left <- read.csv(paste0("/.../left_gbj_pvalue",i,".csv"))
  # gbj.pvalue.list.left[[i]] <- gbj.pvalue.load.left$GBJ.pvalue
  # gbj.pvalue.list.left[[i]][gbj.pvalue.list.left[[i]]==1] <- 0.99
  # p.value = CCT(as.matrix(gbj.pvalue.list.left[[i]]))
  # cct.gbj.pvalue.left[i,] <- c(i,p.value)
  
  gbj.pvalue.load.right <- read.csv(paste0("/.../right_gbj_pvalue",i,".csv"))
  gbj.pvalue.list.right[[i]] <- gbj.pvalue.load.right$GBJ.pvalue
  gbj.pvalue.list.right[[i]][gbj.pvalue.list.right[[i]]==1] <- 0.99
  p.value = CCT(as.matrix(gbj.pvalue.list.right[[i]]))
  cct.gbj.pvalue.right[i,] <- c(i,p.value)
  print(c(i,cct.gbj.pvalue.left[i,2],cct.gbj.pvalue.right[i,2]))
}
# write.csv(cct.gbj.pvalue.left, paste0("/.../adamant_global_left.csv"), row.names = F)
write.csv(cct.gbj.pvalue.right, paste0("/.../adamant_global_right.csv"), row.names = F)



