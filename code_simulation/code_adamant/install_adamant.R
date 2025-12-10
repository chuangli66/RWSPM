install.packages("adamant") # package ‘adamant’ is not available (for R version 3.6.1)
install.packages("devtools")
library(devtools)
install_github("dspluta/adamant")
library(adamant)

# example
X <- matrix(rnorm(100), nrow = 20, ncol = 5)
Y <- rnorm(20)
lambdas_X <- c(0, 1, 1000, Inf)
res = adamant(X, Y, lambdas_X)

# check the paper what is lambda