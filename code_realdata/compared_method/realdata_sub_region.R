rm(list=ls())

width = 3
idx1 = slide_width_step(width/2,width)
idx = idx1$idx
load("./data/left_res_ori.Rdata")

t = seq(0,1,0.05)

m = nrow(idx)
n = nrow(left_res_hipp)

sub_region = list()
for(j in 1:m){
  print(j)
   sub_y = matrix(NA, nrow = n, ncol = 9)
  for(i in 1:n){
    sub_y[i,] = left_res_hipp[i,idx[j,]]
  }
  sub_region[[j]] = sub_y
}

save(sub_region, file = "./data/left_res_hipp_width3.Rdata")


