library('dplyr')


est_drig <- function(data, gamma, y_idx=1, del_idx=NULL, unif_weight=FALSE) {
  if (is.null(del_idx)) {
    del_idx <- y_idx
  }
  m <- length(data)
  if (unif_weight) {
    w <- rep(1/m, m)
  } else {
    w <- sapply(data, function(e) nrow(e))
    w <- w / sum(w)
  }
  
  gram_x <- list()
  gram_xy <- list()
  
  for (e in 1:m) {
    data_e <- data[[e]]
    n <- nrow(data_e)
    y <- data_e[, y_idx]
    x <- data_e[, -c(y_idx, del_idx)]
    
    gram_x[[e]] <- t(x) %*% x / n
    gram_xy[[e]] <- t(x) %*% y / n
  }
  G <- (1 - gamma) * gram_x[[1]] + gamma * Reduce('+', Map('*', gram_x, w))
  
  
  
  Z <- (1 - gamma) * gram_xy[[1]] + gamma * Reduce('+', Map('*', gram_xy, w))
  drig= solve(G) %*% Z

  if(is.infinite(gamma)) {
    G <-  gram_x[[1]] - gram_x[[2]]
    Z <- gram_xy[[1]] - gram_xy[[2]]
    dantzig=solve(G) %*% Z
    return(dantzig)
  }
  return(drig)
}
