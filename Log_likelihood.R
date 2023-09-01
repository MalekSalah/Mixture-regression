Log_lik_p <- function (p, Theta1, Theta2, Y){
  res <- vector(mode = "numeric", length(p))
  for (i in 1:length(p)){res[i] <- sum(log(p[i] * dnorm(Y, Theta1$mu, Theta1$sigma) + (1-p[i]) * dnorm(Y, Theta2$mu, Theta2$sigma)))}
  return(res)
}

Log_lik_etaPhi_p <- function (p, etaPhi, sigma, mu,  Y){
  mu1 <- mu - sigma * etaPhi[3] * (sqrt(1-p)/sqrt(p))
  mu2 <- mu + sigma * etaPhi[3] * (sqrt(p)/sqrt(1-p))
  sigma1 <- sigma * (etaPhi[1] / sqrt(p))
  sigma2 <- sigma * (etaPhi[2] / sqrt(1-p))
  res <- sum(log(p*dnorm(Y, mu1, sigma1)  + (1-p)* dnorm(Y, mu2, sigma2)))
  return(res)
}

Log_lik <- function (p, etaPhi, beta, sigma, X, Y){
  # Compute expectation
  mu_XB <- X %*% beta 
  
  # Components' sd and expectation
  mu1 <- mu_XB - sigma * etaPhi[3] * (sqrt(1-p)/sqrt(p))
  mu2 <- mu_XB + sigma * etaPhi[3] * (sqrt(p)/sqrt(1-p))
  sigma1 <- sigma * (etaPhi[1] / sqrt(p))
  sigma2 <- sigma * (etaPhi[2] / sqrt(1-p))
  res <- sum(log(p*dnorm(Y, mu1, sigma1)  + (1-p)* dnorm(Y, mu2, sigma2)))  
  return(res)
}