

quick_metafor <- function(yi, sei){
  yi <- yi[is.finite(sei)]
  sei <- sei[is.finite(sei)]


  n <- length(yi)
  if(n < 2){
    return(list(b = NA_real_, tau2 = NA_real_, se = NA_real_, pval = NA_real_))
  }
  ytilde <- yi - mean(yi)
  tau2 <- max(0, (sum(ytilde^2) - sum((1 - 1/n) * sei^2)) / (n-1))
  old_tau2 <- tau2
  iter <- 1
  converged <- FALSE
  while(iter <= 10 && ! converged){
    wi   <- 1/(sei^2 + tau2)
    sum_w <- sum(wi)
    sum_w2 <- sum(wi^2)
    sum_w3 <- sum(wi^3)
    weighted_mean <- sum(ytilde * wi) / sum_w
    weighted_RSS <- sum((wi * (ytilde - weighted_mean))^2)
    adj <- (weighted_RSS - sum_w + sum_w2 / sum_w) / (sum_w2 - 2 * sum_w3 / sum_w + (sum_w2 / sum_w)^2)
    tau2 <- max(0,tau2 + adj)
    converged <- abs(1 - (tau2 + 0.5) / (old_tau2 + 0.5)) < 1e-6
    old_tau2 <- tau2
    iter <- iter + 1
  }
  sum_w <- sum(wi)
  weighted_mean <- sum(ytilde * wi) / sum_w
  b <- mean(yi) + weighted_mean

  se <- sqrt(1/sum_w)
  pval <- 2 * min(pnorm(b / se, lower.tail = FALSE),
                  pnorm(b / se, lower.tail = TRUE))
  list(b = b, tau2 = tau2, se = se, pval = pval)
}
