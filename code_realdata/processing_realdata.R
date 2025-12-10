rm(list = ls()); 
wd="" ########please set working directory as needed to store the data to be used in later analysis.
setwd(wd)
source("./code_RWSPM/slide_width.R")
load(file = "./data/left_hipp_ori.Rdata")
load(file = "./data/right_hipp_ori.Rdata")
load(file = "./data/covarites_inf.Rdata")

#####################################################################################################################
########## Remove the effects of baseline covariates on the LQD values of each subregion in the original image ######
#####################################################################################################################
slide_width_value = 11
left_idx_list = slide_width_step(b=ceiling(slide_width_value/4),width = slide_width_value)
leftG = LQD(left_hipp_inf1_3,left_idx_list)
right_idx_list = slide_width_step(b=ceiling(slide_width_value/4),width = slide_width_value)
rightG = LQD(right_hipp_inf1_3,right_idx_list)

XX <- cov_inf1_3%*%solve(t(cov_inf1_3)%*%cov_inf1_3)%*%t(cov_inf1_3)
right_res_G = list()
left_res_G = list()
for (i in seq(length(rightG))) {
  right_res_G[[i]] <- as.matrix(rightG[[i]]) - as.matrix(XX)%*%as.matrix(rightG[[i]])
  print(i)
}
for (i in seq(length(leftG))) {
  left_res_G[[i]] <- as.matrix(leftG[[i]]) - as.matrix(XX)%*%as.matrix(leftG[[i]])
  print(i)
}

save(right_res_G,file = "./data/right_res_G.Rdata")
save(left_res_G,file = "./data/left_res_G.Rdata")


#############################################################################
########## remove the effects of baseline covariates on original image ######
#############################################################################
XX <- cov_inf1_3%*%solve(t(cov_inf1_3)%*%cov_inf1_3)%*%t(cov_inf1_3)
left_res_ori <- as.matrix(left_hipp_inf1_3) - as.matrix(XX)%*%as.matrix(left_hipp_inf1_3)
right_res_ori <- as.matrix(right_hipp_inf1_3) - as.matrix(XX)%*%as.matrix(right_hipp_inf1_3)

save(left_res_ori, file = "./data/left_res_ori.Rdata")
save(right_res_ori, file = "./data/right_res_ori.Rdata")


#############################################################################
########## save left_res_ori and  right_res_ori as .mat file ################
#############################################################################
## install.packages("R.matlab")
library(R.matlab)
writeMat("./data/left_res_ori.mat",  left_res_ori  = as.matrix(left_res_ori))
writeMat("./data/right_res_ori.mat", right_res_ori = as.matrix(right_res_ori))


#############################################################################
########## load the SNPs information (chr SNPs position Allele1 Allele2) ####
#############################################################################
#load("./data/snp_inf1_3.RData")


#############################################################################
########## save covarites_inf.Rdata as covariates_inf.csv ####
#############################################################################
# load(file = "./data/covarites_inf.Rdata")
# colnames(cov_inf1_3) <- c("GENDER","AGE","pc1","pc2","pc3","pc4","pc5","APOE-4")
# write.csv(cov_inf1_3, file = "./data/covarites_inf.csv")
