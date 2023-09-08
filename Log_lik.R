Log_lik <- function (p_val, etaPhi_val, sigma_val, mu_val, dt){
  mu1 <- mu_val - sigma_val * etaPhi_val[3] * (sqrt(1-p_val)/sqrt(p_val))
  mu2 <- mu_val + sigma_val * etaPhi_val[3] * (sqrt(p_val)/sqrt(1-p_val))
  sigma1 <- sigma_val * (etaPhi_val[1] / sqrt(p_val))
  sigma2 <- sigma_val * (etaPhi_val[2] / sqrt(1-p_val))
  res <- sum(log(p_val*dnorm(dt, mu1, sigma1)  + (1-p_val)* dnorm(dt, mu2, sigma2)))
  return(res)
}
