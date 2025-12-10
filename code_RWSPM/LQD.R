LQD <- function(y, idx, t = seq(0,1,0.05)){

  m = nrow(idx) 
  Y = as.matrix(y); rm(y)
  n = nrow(Y) 
  
  
  QD = list()
  for(j in 1:m){
    # print(j)
    qd = matrix(NA, nrow = n, ncol = length(t))
    for(i in 1:n){
      y = Y[i,idx[j,]]
      
      d = density(y)
      d$y[d$y==0]=0.0001
      
      # QD
      q = quantile(y, t)
      qd[i,] = approx(d$x, d$y, q)$y
    }
    QD[[j]] = qd
  }
  
  G = list()
  for(j in 1:m){
    G[[j]] = -log(QD[[j]])
  }
  return(G)
}


