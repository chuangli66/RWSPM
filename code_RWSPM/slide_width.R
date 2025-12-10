library(overlapping)
# Calculate the JSD between two kernel density estimation functions
jsd_from_kde_samples <- function(dens_x, densy) {
  # KDE 密度估计
  dx <- dens_x$x[2] - dens_x$x[1]

  # 转为概率分布（归一化）
  p <- dens_x$y / sum(dens_x$y * dx)
  q <- dens_y$y / sum(dens_y$y * dx)

  # 加平滑项，避免 log(0)
  eps <- 1e-10
  p <- p + eps
  q <- q + eps
  p <- p / sum(p)
  q <- q / sum(q)

  # 平均分布
  m <- 0.5 * (p + q)

  # KL 散度部分
  kl <- function(a, b) sum(a * log(a / b))

  jsd <- 0.5 * kl(p, m) + 0.5 * kl(q, m)
  return(jsd)
}
# The golden section method is used to find the optimal slide width that minimizes mean(ov2).
golden_section_search <- function(y, f, a, b, M1, M2, n_sam, tol = 1) {
  phi <- (1 + sqrt(5)) / 2
  resphi <- 1/ phi

  x1 <- b - (b - a) * resphi
  x2 <- a + (b - a) * resphi
  f1 <- f(y, x1, M1, M2)
  f2 <- f(y, x2, M1, M2)

  while (abs(x1 - x2) > tol) {
    if (f1 < f2) {
      b <- x2
      x2 <- x1
      f2 <- f1
      x1 <- b - (b - a) * resphi
      f1 <- f(y, x1, M1, M2)
    } else {
      a <- x1
      x1 <- x2
      f1 <- f2
      x2 <- a + (b - a) * resphi
      f2 <- f(y, x2, M1, M2)
    }
    print(c(abs(x1 - x2),1))
  }
  return(ceiling((a + b) / 2))
}

# Objective function: Given a sliding window width w, calculate mean(ov2).
compute_ov <- function(y, window_width, M1, M2, n_grid = 256) {
  n_sam <- nrow(y)
  window_width <- ceiling(window_width)
  w1 <- ceiling(window_width / 4)
  idx1 <- slide_width_step(w1, window_width, M1, M2)

  ov2 <- numeric(nrow(idx1$idx))

  # The KDE function
  kde_prob <- function(x, grid) {
    dens <- density(x, from = min(grid), to = max(grid), n = length(grid))
    p <- dens$y / sum(dens$y)
    return(p)
  }

  for (h in seq_len(nrow(idx1$idx))) {
    sub_region_idx <- idx1$idx[h, ]
    sub_region <- y[, sub_region_idx, drop = FALSE]  # n_sam x window_width

    # Unified Grid
    grid <- seq(min(sub_region), max(sub_region), length.out = n_grid)

    # KDE matrix
    prob_mat <- t(apply(sub_region, 1, kde_prob, grid = grid))  # n_sam x n_grid

    # JSD matrix
    # 1. Construct a broadcast matrix for all rows.
    P <- array(rep(prob_mat, each = n_sam), dim = c(n_sam, n_sam, n_grid))
    Q <- aperm(P, c(2, 1, 3))
    M <- 0.5 * (P + Q)

    # 2. KL(P||M) and KL(Q||M)
    eps <- 1e-10
    P <- P + eps
    Q <- Q + eps
    M <- M + eps

    KL_PM <- rowSums(P * log(P / M), dims = 2)  # n_sam x n_sam
    KL_QM <- rowSums(Q * log(Q / M), dims = 2)

    jsd_matrix <- 0.5 * (KL_PM + KL_QM)

    # ov2
    ov2[h] <- mean(jsd_matrix)
    # print(c(h,ov2[h]))
  }

  ov_mean <- mean(ov2)
  ov_loss <- ov_mean + 0.6 * (window_width - 2) / (ceiling(min(M1,M2)/2) - 2)
  print(paste0("width_", window_width, ": ", ov_mean,"-",ov_loss))
  return(ov_loss)
}

# Based on the set sliding step size b and sliding window width, segment the image of size M1*M2.
slide_width_step <- function(b, width, M1 = 150, M2 = 100){
  g0 = matrix(NA, nrow = M1, ncol = M2)
  idx = matrix(NA, ncol = width*width, nrow = (ceiling(M1/b) - 1)*(ceiling(M2/b) - 1))

  x1 = 1            # start point in x
  x2 = x1 + width - 1 # end point in x

  i = 1
  i1 = 0
  final_x = 0
  width_x = 1
  width_y = 1
  while(TRUE){
    y1 = 1            # start point in y
    y2 = y1 + width - 1 # end point in y
    while(TRUE){
      if(y2 <= M2){
        # find the idx in 1:15000
        g = g0
        g[x1:x2, y1:y2] = T
        gl = as.vector(t(g)) # default: stack all the columns
        idx[i,] = which(gl)
        # idx = rbind(idx,which(gl))
        y1 = y1 + b
        y2 = y1 + width - 1
        i = i+1
        # print(y2)
        if(y2 == M2 + b){
          break
        }
      }else{
        y2 = M2
        y1 = y2 - width + 1 # 控制行扫描超出图像边界时候的子区域
        g = g0
        g[x1:x2, y1:y2] = T
        gl = as.vector(t(g))
        idx[i,] = which(gl)
        i = i + 1
        if(y2 == M2 & x2 == M1){
          return(list("width_x" = width_x,
                      "width_y" = dim(na.omit(idx))[1]/width_x,
                      "idx" = na.omit(idx)))
        }
        break
      }
    }

    i1 = i1+1
    if(x2 == M1 || final_x == 1){
      return(list("width_x" = width_x,
                  "width_y" = dim(na.omit(idx))[1]/width_x,
                  "idx" = na.omit(idx)))
    }
    x1 = x1 + b
    x2 = x1 + width - 1
    width_x = width_x + 1
    # print(c(x1,x2))

    if(x2 > M1){
      x2 = M1
      x1 = M1 - width + 1
      final_x = 1
    }
  }
}


