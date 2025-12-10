# reformat G into a 3D matrix and save (as R object)
# library(abind)
# install.packages("RcppCNPy")
library(RcppCNPy)

## k1 is the number of Setting 1-7
k1=1
## File loading path
root = paste0("/.../Setting",k1,"/")
for(l in 1:100){
  print(l)
  p = paste0(root,l,'/')
  fn = paste0(p,"Y_slide", l, ".RData")
  print(fn)
  if(file.exists(fn)){
    load(fn)
    Y_slide = do.call(rbind, Y_slide) 
    # save(LQD, file=paste0(p,"LQD_3Dmat.RData")) # can't be read into python using pyreadr
    npySave(paste0(p,"Y_slide.npy"), Y_slide) # doesn't work for 3d matrix
  } else {
    print("file doesn't exist!")
  }
}

