## instal package RWSPM (/.../RWSPM_0.1.0.tar.gz)
library(RWSPM)
## File loading path
root = "/.../Setting" 

## regional partition, Window width set to 3, step size set to 1.5.
idx_list = slide_width_step(b = 1.5, width = 3)
subr_idx = idx_list$idx
for (k in seq(7)) {
  # k=4
  cat(paste0("set",k),"\n")
  for (l in seq(1,100)) {
    cat(l,"\n")
    #
    load(paste0(root, k, "/sim_", l , "/y.RData"))
    m = nrow(subr_idx) # number of (sub)regions: 126
    Y = as.matrix(y); rm(y)
    n = nrow(Y) # number of subjects (sample size): 200

    # file_path = paste0(root, k, "sim_", l ,"/Y_slide",l,".RData")

    Y_slide = list()
    for(j in 1:m){
      # print(j)
      qd = matrix(NA, nrow = n, ncol = dim(subr_idx)[2])
      for(i in 1:n){
        qd[i,] = Y[i,subr_idx[j,]]
      }
      Y_slide[[j]] = qd
    }
    save(Y_slide, file = paste0(root, k, "/sim_", l ,"/Y_slide",l,".RData"))
  }
}

library(R.matlab)
### Load the segmented images as a .csv file for running the fvgwas method.
for (k in seq(7)) {
  root = "./set"
  cat(paste0("set",k),"\n")
  for (l in seq(100)){
    cat(l,"\n")
    load(paste0(root, k, "/sim_", l , "/y.RData"))
    m = nrow(subr_idx) # number of (sub)regions: 126
    Y = as.matrix(y); rm(y)
    n = nrow(Y) # number of subjects (sample size): 200

    Y_slide_matrix = matrix(NA,dim(subr_idx)[1],200*dim(subr_idx)[2])

    folder_path = paste0(root, k, "/sim_", l ,"/Y_slide",l,".mat")
    if (!dir.exists(folder_path)) {

      Y_slide = list()
      for(j in 1:m){
        # print(j)
        qd = matrix(NA, nrow = n, ncol = dim(subr_idx)[2])
        for(i in 1:n){
          qd[i,] = Y[i,subr_idx[j,]]
        }
        Y_slide[[j]] = qd
      }
      for (m in seq(length(Y_slide))) {
        Y_slide_matrix[m,] = as.vector(Y_slide[[m]])
      }
      writeMat(paste0(root, k, "/sim_", l ,"/Y_slide",l,".mat"), Y = Y_slide_matrix)
      cat("The file has been created.","\n")
    } else {
      cat("The file already exists.:", folder_path, "\n")
    }
  }
}



