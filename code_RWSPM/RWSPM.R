source("./code_RWSPM/LQD.R")
RWSPM <- function(x,
                  y,
                  M1,
                  M2,
                  window.width = NA,
                  method = "bcov"){
  ##################################
  # y: High-dimensional image data
  # x: Predictive variables (e.g. gene data)
  # M1: The height of the image
  # M2: The width of the image
  # window.width: Then width of sliding window
  # method: The selected test statistic ("bcov", "dcov")
  ##################################
  # Selecting Window Width Using the Golden Section Method
  y = as.matrix(y)
  dim_y = dim(y)
  if (is.na(window.width)){
    cat("Regional Partitioning")
    optimal_width = golden_section_search(y, compute_ov, a = 2, b = ceiling(min(M1,M2)/2), M1 = M1, M2 = M2, n_sam = dim_y[1])
    # cat("The optimal width of the sliding window is：", round(optimal_width), "\n")
  }else{
    optimal_width = window.width
  }
  # The partition of the image
  idx1 = slide_width_step(b = ceiling(window_width / 4), width = optimal_width, M1 = M1, M2= M2)
  cat("LQD transformation for each sub-region\n")
  # LQD transformation for each sub-region
  G = LQD(y,idx1$idx)

  m = nrow(idx1$idx) # The number of sub-region
  n = nrow(y) # sample size

  rw.pvalue = matrix(NA, nrow = m, ncol = 3)
  colnames(rw.pvalue) = c("region.j",
                       method,
                       "pvalue")
  cat("Running region-wise independence test across all subregions\n")
  if(method == "bcov"){
    library(Ball)
    for(j in 1:m){
      res = bcov.test(x = x, y = G[[j]], method = "limit", weight = "constant")
      rw.pvalue[j,] = c(j,
                          res$statistic,
                          res$p.value)
    }
  }
  if(method == "dcov"){
    library(dcortools)
    for(j in 1:m){
      res = distcov.test(as.matrix(x),as.matrix(G[[j]]), method = "gamma")
      rw.pvalue[j,] = c(j,
                          res$dcov,
                          res$pvalue)
    }
  }
  return(list(rw.pvalue = rw.pvalue,
              num.region.r = idx1$width_y,
              num.region.c = idx1$width_x))
}

CCT_chisq <- function(pvals,weights=NULL){
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

RW.MTCCT <- function(pvalue_mat, num.region.r, num.region.c, threshold = 0.05){
  ## num.region.r corresponds to M1 of the image (M2 = 100)
  ## num.region.c corresponds to M2 of the image (M1 = 150)
  # num.region.r = idx$num.region.r
  # num.region.c = idx$num.region.c
  library(igraph)
  rw.mtcct <- NULL
  if(is.null(dim(pvalue_mat))){
    pvalue_mat = as.matrix(pvalue_mat)
  }

  for (boot.i in seq(dim(pvalue_mat)[1])) {
    # p-value of the boot.i row
    pvec <- pvalue_mat[,boot.i]

    # 1. the position of the smallest significant p-value
    min_sig <- min(pvec[pvec < threshold])
    min_pos <- which(pvec == min_sig)

    # Convert to matrix coordinates (row, col)
    min_x <- ((min_pos - 1) %/% num.region.r) + 1
    min_y <- ((min_pos - 1) %%  num.region.r) + 1

    # 2. Construct a significant matrix
    o.matrix <- t(matrix(pvec, num.region.r, num.region.c))
    o.index <- which(o.matrix < threshold, arr.ind = TRUE)

    # Calculate the distance matrix between the significant points
    o.dist <- as.matrix(dist(o.index))

    # Construct a graph (adjacency relationship: distance of 1)
    graph <- graph_from_adjacency_matrix(o.dist == 1, mode = "undirected")

    # Search for connected regions
    connected_components <- clusters(graph)


    # Find the component number containing min_pos
    # First, locate the positions of min_x and min_y in o.index
    selected_idx <- which(o.index[,1] == min_x & o.index[,2] == min_y)

    comp_id <- connected_components$membership[selected_idx]

    # Extract the indices (row numbers in o.index) of all the points in this connected region
    comp_vertices <- which(connected_components$membership == comp_id)

    # 3. Restore the original index from the coordinates: detecte_index
    detecte_index <- (o.index[comp_vertices,1] - 1) * num.region.r + o.index[comp_vertices,2]

    # Remove the sub-significant area index
    non_sig_count <- sum(pvec >= threshold)
    cct_region_count = non_sig_count + length(detecte_index)

    # Obtain the CCT p-value of this connected region
    if(length(detecte_index) == 0){return(1)}else{
      weight_cct <- rep(0,length(pvec))
      weight_cct[detecte_index] <- 1/cct_region_count

      rw.mtcct <- CCT_chisq(pvec,weights = weight_cct)
      # CCT_chisq(pvec,weights = weight_cct)
      # CCT(pvec,weights = weight_cct)

      return(rw.mtcct)}
  }
}


