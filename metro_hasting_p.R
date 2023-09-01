Log_post_p <- function (p, Theta1, Theta2, Tau, Y) {
  return(Log_lik(p, Theta1, Theta2, Y) + log(dbeta(p, Tau$tau1, Tau$tau2)))
}

# This is the vectorized version of the Metropolis Hastings algorithm used to sample p
# P_0, the starting point could be a vector 
metro_hasting_p <- function(P_0, eps_p, T){
  P <- matrix(NA, nrow = T, ncol = length(P_0))
  P[1,] <- P_0
  acceptance_bool <- matrix(nrow = T-1, ncol = length(P_0))
  for (i in 1:(T-1)){
    # Generate proposal
    p_t <- P[i,]
    p_prop <- rbeta(length(p_t), p_t*eps_p+1, (1-p_t)*eps_p+1)
    
    #computing acceptance probability
    log_acceptance <- (Log_post_p(p_prop, Theta1, Theta2, Tau, Y)+log(dbeta(p_t, p_prop*eps_p+1, (1-p_prop)*eps_p+1))) - (Log_post_p(p_t, Theta1, Theta2, Tau, Y) + log(dbeta(p_prop, p_t*eps_p+1, (1-p_t)*eps_p+1)))
    acc <- vector(mode = "numeric", length = length(log_acceptance))
    for (k in 1:length(log_acceptance)) acc[k] <- min(1, exp(log_acceptance[k]))
    
    #Decision on the next step
    U <- runif(length(acc))
    for (k in 1:length(acc)){
      if ( acc[k] > U[k]) {
        P[i+1,k] <- p_prop[k]
      }
      else P[i+1,k] <- p_t[k]
      acceptance_bool[i,] <- (acc > U)
    }
  }
#Acceptance rate
acceptance_rate <- colSums(acceptance_bool)/T
return(list("acceptance_rate" = acceptance_rate, "Sample" = P))
}
