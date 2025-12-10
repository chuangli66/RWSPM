## instal package RWSPM (/.../RWSPM_0.1.0.tar.gz)
library(RWSPM)
library(igraph)
library(ggplot2)

idx1 = slide_width_step(ceiling(10/4),10)
significant_overlapping_region_indices_width10 = idx1$idx

idx1 = slide_width_step(1.5,3)
significant_overlapping_region_indices_width3 = idx1$idx


cover_function_bcov1 <- function(bcov.pvalue1){
  cover_size <- NULL
  bcov.pvalue1.bh <- matrix(NA,dim(bcov.pvalue1)[1],dim(bcov.pvalue1)[2])
  for (boot.i in seq(dim(bcov.pvalue1)[1])) {
    o.bcov.matrix <- t(matrix(bcov.pvalue1[boot.i,],19,29))
    o.bcov.index <- which(o.bcov.matrix < 0.05, arr.ind = TRUE)
    inter_size <- intersect(significant_overlapping_region_indices_width10,(o.bcov.index[,1] - 1)*19 + o.bcov.index[,2])
    cover_size[boot.i] <- 2*length(inter_size)/(length(significant_overlapping_region_indices_width10) + length((o.bcov.index[,1] - 1)*19 + o.bcov.index[,2]))
    cover_size[boot.i]
  }
  return(cover_size)
}
un_cover_function <- function(bcov.pvalue1){
  cover_size <- NULL
  bcov.pvalue.bh <- matrix(NA,dim(bcov.pvalue1)[1],dim(bcov.pvalue1)[2])
  for (boot.i in seq(dim(bcov.pvalue)[1])) {
    index1 <- which(bcov.pvalue1[boot.i,] < 0.05, arr.ind = TRUE)
    inter_size <- intersect(significant_overlapping_region_indices_width10,index1[,2])
    union_size <- union(significant_overlapping_region_indices_width10,index1[,2])
    FP = (length(index1[,2]) - length(inter_size))
    TN = 551 - length(union_size)
    cover_size[boot.i] <- FP/(TN+FP)
  }
  return(cover_size)
} 

cover_function_compare1 <- function(bcov.pvalue1){
  cover_size <- NULL
  bcov.pvalue1.bh <- matrix(NA,dim(bcov.pvalue1)[1],dim(bcov.pvalue1)[2])
  for (boot.i in seq(dim(bcov.pvalue1)[1])) {
    o.bcov.matrix <- t(matrix(bcov.pvalue1[boot.i,],66,99))
    o.bcov.index <- which(o.bcov.matrix < 0.05, arr.ind = TRUE)
    
    inter_size <- intersect(significant_overlapping_region_indices_width3,(o.bcov.index[,1] - 1)*66 + o.bcov.index[,2])
    cover_size[boot.i] <- 2*length(inter_size)/(length(significant_overlapping_region_indices_width3) + length((o.bcov.index[,1] - 1)*66 + o.bcov.index[,2]))
    cover_size[boot.i]
  }
  return(cover_size)
}
un_cover_function_compare <- function(bcov.pvalue1){
  cover_size <- NULL
  bcov.pvalue.bh <- matrix(NA,dim(bcov.pvalue1)[1],dim(bcov.pvalue1)[2])
  for (boot.i in seq(dim(bcov.pvalue1)[1])) {
    index1 <- which(bcov.pvalue1[boot.i,] < 0.05, arr.ind = TRUE)
    inter_size <- intersect(significant_overlapping_region_indices_width3,index1[,2])
    union_size <- union(significant_overlapping_region_indices_width3,index1[,2])
    FP = (length(index1[,2]) - length(inter_size))
    TN = 6534 - length(union_size)
    cover_size[boot.i] <- FP/(TN+FP)
  }
  return(cover_size)
} 

cover_bcov <- data.frame()
cover_gbj <- data.frame()
cover_adamant <- data.frame()
cover_rdcov <- data.frame()
cover_skpcr <- data.frame()
cover_fvgwas <- data.frame()
uncover_bcov <- data.frame()
uncover_gbj <- data.frame()
uncover_adamant <- data.frame()
uncover_rdcov <- data.frame()
uncover_skpcr <- data.frame()
uncover_fvgwas <- data.frame()

## k1 is changed from 2 to 7
k1 = 2
save_root = paste0("./.../Setting",k1,"/")
### bcov
bcov.pvalue <- read.table(file = paste0(save_root,"pval_bcov_set",k1,".csv") ,sep = ",",header=TRUE)
### gbj
gbj.pvalue <- read.table(file = paste0(save_root,"pval_gbj_set",k1,".csv"),sep = ",",header=TRUE)
### adamant
adamant.pvalue <- read.table(file = paste0(save_root,"pval_adamant_set",k1,".csv"),sep = ",",header=TRUE)
### rdcov
rdcov.pvalue <- read.table(file = paste0(save_root,"pval_rdcov_set",k1,".csv"),sep = ",",header=TRUE)

### skpcr 
skpcr.pvalue <- read.table(file = paste0(save_root,"pval_skpcr_set",k1,".csv"),sep = ",",header=TRUE)
### fvgwas
fvgwas.pvalue <- read.table(file = paste0(save_root,"_pval_fvgwas_set",k1,".csv"),sep = ",",header = TRUE)



### 
uncover.bcov.pvalue1 <- un_cover_function(bcov.pvalue[,-1])
uncover.gbj.pvalue1 <- un_cover_function_compare(gbj.pvalue[,-1])
uncover.adamant.pvalue1 <- un_cover_function_compare(adamant.pvalue[,-1])
uncover.rdcov.pvalue1 <- un_cover_function_compare(rdcov.pvalue[,-1])
uncover.skpcr.pvalue1 <- un_cover_function_compare(skpcr.pvalue[,-1])
uncover.fvgwas.pvalue1 <- un_cover_function_compare(fvgwas.pvalue[,-1])


cover.bcov.pvalue1 <- cover_function_bcov1(bcov.pvalue[,-1])
cover.gbj.pvalue1 <- cover_function_compare1(gbj.pvalue[,-1])
cover.adamant.pvalue1 <- cover_function_compare1(adamant.pvalue[,-1])
cover.rdcov.pvalue1 <- cover_function_compare1(rdcov.pvalue[,-1])
cover.skpcr.pvalue1 <- cover_function_compare1(skpcr.pvalue[,-1])
cover.fvgwas.pvalue1 <- cover_function_compare1(fvgwas.pvalue[,-1])

print(cover.bcov.pvalue1)
print(cover.gbj.pvalue1)
print(cover.adamant.pvalue1)
print(cover.rdcov.pvalue1)
print(cover.skpcr.pvalue1)
print(cover.fvgwas.pvalue1)


cover_mat <- c(cover.bcov.pvalue1,cover.gbj.pvalue1,cover.adamant.pvalue1,cover.rdcov.pvalue1,cover.skpcr.pvalue1,cover.fvgwas.pvalue1)
uncover_mat <- c(uncover.bcov.pvalue1,uncover.gbj.pvalue1,uncover.adamant.pvalue1,uncover.rdcov.pvalue1,uncover.skpcr.pvalue1,uncover.fvgwas.pvalue1)

Method  <- c(rep("RW-SPM",length(cover.bcov.pvalue1)),rep("GBJ",length(cover.gbj.pvalue1)),rep("AdaMant",length(cover.adamant.pvalue1)),rep("Rdcov",length(cover.rdcov.pvalue1)),rep("sKPCR",length(cover.skpcr.pvalue1)),rep("fvGWAS",length(cover.fvgwas.pvalue1)))
unMethod  <- c(rep("RW-SPM",length(uncover.bcov.pvalue1)),rep("GBJ",length(uncover.gbj.pvalue1)),rep("AdaMant",length(uncover.adamant.pvalue1)),rep("Rdcov",length(uncover.rdcov.pvalue1)),rep("sKPCR",length(uncover.skpcr.pvalue1)),rep("fvGWAS",length(uncover.fvgwas.pvalue1)))

cover_dataframe <- data.frame(Cover_rate = cover_mat,Method = Method)
uncover_dataframe <- data.frame(unCover_rate = uncover_mat,unMethod = unMethod)

cover_dataframe$Method = factor(cover_dataframe$Method, levels = c("RW-SPM","GBJ","AdaMant","Rdcov","sKPCR","fvGWAS"))
uncover_dataframe$unMethod = factor(uncover_dataframe$unMethod, levels = c("RW-SPM","GBJ","AdaMant","Rdcov","sKPCR","fvGWAS"))


ggplot(cover_dataframe, aes(x = Method, y = cover_mat, fill = Method)) +
  geom_boxplot()+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  # theme(legend.background = element_rect(fill = "lightgrey"))+
  labs(title = "(a)", x = "method", y = "Dice overlap ratio (DOR) ")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )+
  ylim(c(0,1))

ggplot(uncover_dataframe, aes(x = unMethod, y = unCover_rate, fill = unMethod)) +
  geom_boxplot()+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  # theme(legend.background = element_rect(fill = "lightgrey"))+
  labs(title = "(a)", x = "method", y = "False positive rate (FDR)")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )+
  ylim(c(0,1))

