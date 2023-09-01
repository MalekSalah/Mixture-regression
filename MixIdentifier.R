MixIdentifier <- function(Y, T){
  T_adapt <- 300
  T <- T+100
  
  #Log-likelihood function
  Log_lik <- function (p_val, etaPhi_val, sigma_val, mu_val, dt){
    mu1 <- mu_val - sigma_val * etaPhi_val[3] * (sqrt(1-p_val)/sqrt(p_val))
    mu2 <- mu_val + sigma_val * etaPhi_val[3] * (sqrt(p_val)/sqrt(1-p_val))
    sigma1 <- sigma_val * (etaPhi_val[1] / sqrt(p_val))
    sigma2 <- sigma_val * (etaPhi_val[2] / sqrt(1-p_val))
    res <- sum(log(p_val*dnorm(dt, mu1, sigma1)  + (1-p_val)* dnorm(dt, mu2, sigma2)))
    return(res)
  }
  
  #Log-posterior distributions:
  Log_post_p <- function (p_val, etaPhi_val, sigma_val, mu_val, Tau, dt) {
    return( Log_lik(p_val, etaPhi_val, sigma_val, mu_val, dt) + log(dbeta(p_val, Tau$tau1, Tau$tau2)))
  }
  Log_post_etaPhi <- function (p_val, etaPhi_val, sigma_val, mu_val, Alpha, dt){
    return( Log_lik(p_val, etaPhi_val, sigma_val, mu_val, dt) + ddirichlet(matrix(etaPhi_val, nrow = 1), Alpha, log = TRUE))
  }
  Log_post_sigma <- function ( p_val, etaPhi_val, sigma_val, mu_val, ksi, dt){
    return( Log_lik(p_val, etaPhi_val, sigma_val, mu_val, dt) + dinvgamma(sigma_val, ksi$shape, ksi$rate, log = TRUE))
  }
  Log_post_mu <- function (p_val, etaPhi_val, sigma_val, mu_val, psi, dt){
    return( Log_lik(p_val, etaPhi_val, sigma_val, mu_val, dt) + dnorm(mu_val,psi$mean, psi$sd, log = TRUE ))
  }
  
  #Prior distribution parameters:
  Tau <- list("tau1"=10 , "tau2"=10)
  Alpha <- c(.5,.5,.5)
  ksi <- list("shape"=1.5, "rate"=2)
  psi <- list("mean"=1, "sd"=1)
  
  
  
  # Mixture weights p
  eps_p <- 0.4
  P <- numeric(length = T)
  p_0 <- rbeta(1,.5, .5)
  P[1] <- p_0
  P_adapt <- numeric(length = T_adapt)
  
  # Mixture mean 
  eps_mu <- .4
  Mu <- numeric(length = T)
  Mu[1] <- mean(Y)
  Mu_adapt <- numeric(length = T_adapt)
  
  # Mixture std sigma
  eps_sigma_sqrd <- .4
  sigma_sqrd <- numeric(length = T)
  sigma_sqrd[1] <- sd(Y)
  sigma_sqrd_adapt <- numeric(length = T_adapt)
  
  # Component ecpectations mu1/2
  mu1_0 <- 1
  mu2_0 <- 1
  mu1 <- vector(mode = "numeric", length = T)
  mu2 <- vector(mode = "numeric", length = T)
  mu1[1] <- mu1_0
  mu2[1] <- mu2_0
  
  # Component std sigma1/2
  sigma1_0 <- 1
  sigma2_0 <- 1
  sigma1 <- vector(mode = "numeric", length = T)
  sigma2 <- vector(mode = "numeric", length = T)
  sigma1[1] <- sigma1_0
  sigma2[1] <- sigma2_0
  
  # Initialization etaPhi (contient les carres de eta_1, eta_2 et Phi )
  eps_etaPhi_sqrd <- 300
  etaPhi_sqrd_0 <- rdirichlet(1, Alpha)
  etaPhi_sqrd <- matrix(NA, T, 3)
  etaPhi_sqrd[1, ] <- etaPhi_sqrd_0
  etaPhi_sqrd_adapt <- matrix(NA, T_adapt, 3)
  
  # Algorithm vitals
  ac_p <-  logical(length = T-1)
  ac_rate_p <- numeric(length = T-1)
  
  ac_sigma_sqrd <- logical(length = T-1)
  ac_rate_sigma_sqrd <- numeric(length = T-1)
  
  ac_etaPhi_sqrd <-  logical(length = T-1)
  ac_rate_etaPhi_sqrd <- numeric(length = T-1)
  
  ac_mu <- logical(length = T-1)
  ac_rate_mu <- numeric(length = T-1)
  
  # Adaptation
  ac_rate_lower <- .35
  ac_rate_upper <- .30
  adapt_eps_etaPhi_sqrd <- FALSE
  adapt_eps_p <- FALSE
  adapt_eps_mu <- FALSE
  adapt_eps_sigma_sqrd <- FALSE
  cur_adapt <- 300
  ac_rate_span <- 75
  nb_adapt <- NULL
  

  
  for (i in 1:(T - 1)) {
    # Initialize the Gibbs iteration
    p_t <- P[i]
    etaPhi_sqrd_t <- etaPhi_sqrd[i,]
    sigma_sqrd_t <- sigma_sqrd[i]
    mu_t <- Mu[i]
    ##=======================Sample the mixture weight=====================================>
    #Generate a proposal
    p_prop <- rbeta(1, (p_t*eps_p)+1, ((1-p_t)*eps_p)+1)
    
    # Computing acceptance
    log_acceptance <- 
      (Log_post_p(p_prop, sqrt(etaPhi_sqrd_t), sqrt(sigma_sqrd_t), mu_t, Tau, Y) + dbeta(p_t, p_prop*eps_p+1, (1-p_prop)*eps_p+1, log = TRUE)) -
      (Log_post_p(p_t, sqrt(etaPhi_sqrd_t), sqrt(sigma_sqrd_t), mu_t, Tau, Y) + dbeta(p_prop, p_t*eps_p+1, (1-p_t)*eps_p+1, log = TRUE))
    acc_p <- min(1, exp(log_acceptance))
    
    # Decision on the next step 
    U = runif(1)
    if (s <- (acc_p>U)) p_next <- p_prop
    else p_next <- p_t
    # Acceptance rate
    if (i==1) ac_p[i] <- s
    else ac_p[i] <- ac_p[i-1]+s
    ac_rate_p[i] <- ac_p[i]/i
    
    
    if ((i%%cur_adapt==0) && ((median(ac_rate_p[(i-cur_adapt/2):i])>.4) || (median(ac_rate_p[(i-cur_adapt/2):i])<.3))) {
      i_start <- i
      adapt_eps_p <- TRUE
      ac_rate_window <- numeric(ac_rate_span)
    }
    
    if (adapt_eps_p == TRUE){
      P_adapt[1] <- p_next
      for (j in 1:T_adapt){
        p_t <- P_adapt[j]
        p_prop <- rbeta(1, p_t*eps_p+1, (1-p_t)*eps_p+1)
        
        # Computing acceptance probability
        log_acceptance <- 
          (Log_post_p(p_prop, sqrt(etaPhi_sqrd_t), sqrt(sigma_sqrd_t), mu_t, Tau, Y) + dbeta(p_t, p_prop*eps_p+1, (1-p_prop)*eps_p+1, log = TRUE)) -
          (Log_post_p(p_t, sqrt(etaPhi_sqrd_t), sqrt(sigma_sqrd_t), mu_t, Tau, Y) + dbeta(p_prop, p_t*eps_p+1, (1-p_t)*eps_p+1, log = TRUE))
        acc_p <- min(1, exp(log_acceptance))
        
        # Decision on the next step
        U = runif(1)
        if (s <- (acc_p>U)) P_adapt[j+1] <- p_prop
        else P_adapt[j+1] <- p_t
        
        ac_rate_window[j %% ac_rate_span + 1] <- s
        
        if (j >= ac_rate_span){
          current_ac_rate <- sum(ac_rate_window) / ac_rate_span
          
          # Adjust the eps_p 
          if (current_ac_rate < ac_rate_lower) eps_p <- eps_p*1.1
          else if (current_ac_rate > ac_rate_upper) eps_p <- eps_p*0.9
        }
      }
      adapt_eps_p <- FALSE
    }    
    
    
    P[i+1] <- p_next
    ##======================Sample the eta1, et2 and phi all squared========================>
    # Generate a proposal
    etaPhi_sqrd_prop <-
      rdirichlet(1,
                 c(
                   etaPhi_sqrd_t[1] * eps_etaPhi_sqrd + 1,
                   etaPhi_sqrd_t[2] * eps_etaPhi_sqrd + 1,
                   etaPhi_sqrd_t[3] * eps_etaPhi_sqrd + 1
                 ))
    
    # Computing acceptance probability
    tmp1 <-
      c(
        etaPhi_sqrd_prop[1] * eps_etaPhi_sqrd + 1,
        etaPhi_sqrd_prop[2] * eps_etaPhi_sqrd + 1,
        etaPhi_sqrd_prop[3] * eps_etaPhi_sqrd + 1
      )
    tmp2 <-
      c(etaPhi_sqrd_t[1] * eps_etaPhi_sqrd + 1,
        etaPhi_sqrd_t[2] * eps_etaPhi_sqrd + 1,
        etaPhi_sqrd_t[3] * eps_etaPhi_sqrd + 1)
    
    log_acceptance <-
      (
        Log_post_etaPhi(p_next, sqrt(etaPhi_sqrd_prop), sqrt(sigma_sqrd_t), mu_t, Alpha, Y) + ddirichlet(matrix(etaPhi_sqrd_t, nrow = 1), tmp1 , log = TRUE)
      ) - (Log_post_etaPhi(p_next, sqrt(etaPhi_sqrd_t), sqrt(sigma_sqrd_t), mu_t, Alpha, Y) + ddirichlet(matrix(etaPhi_sqrd_prop, nrow = 1), tmp2, log = TRUE))
    acc_etaPhi_sqrd <- min(1, exp(log_acceptance))
    
    # Decision on the next step
    U = runif(1)
    if (s <- (acc_etaPhi_sqrd > U))
      etaPhi_sqrd_next <- etaPhi_sqrd_prop
    else
      etaPhi_sqrd_next <- etaPhi_sqrd_t
    
    # Acceptance rate
    if (i==1) ac_etaPhi_sqrd[i] <- s
    else ac_etaPhi_sqrd[i] <- ac_etaPhi_sqrd[i-1]+s
    ac_rate_etaPhi_sqrd[i] <- ac_etaPhi_sqrd[i]/i
    
    
    if ((i%%cur_adapt==0) && ((mean(ac_rate_etaPhi_sqrd[(i-cur_adapt/2):i])>.4) || (mean(ac_rate_etaPhi_sqrd[(i-cur_adapt/2):i])<.3))) {
      adapt_eps_etaPhi_sqrd <- TRUE
      nb_adapt <- c(nb_adapt, i)
      ac_rate_window <- numeric(ac_rate_span)
    }
    
    if (adapt_eps_etaPhi_sqrd == TRUE){
      etaPhi_sqrd_adapt[1,] <- etaPhi_sqrd_next
      for (j in 1:(T_adapt-1)){
        etaPhi_sqrd_t <- etaPhi_sqrd_adapt[j,]
        etaPhi_sqrd_prop <-
          rdirichlet(1,
                     c(
                       etaPhi_sqrd_t[1] * eps_etaPhi_sqrd + 1,
                       etaPhi_sqrd_t[2] * eps_etaPhi_sqrd + 1,
                       etaPhi_sqrd_t[3] * eps_etaPhi_sqrd + 1
                     ))
        # Computing acceptance probability
        C1 <-
          c(
            etaPhi_sqrd_prop[1] * eps_etaPhi_sqrd + 1,
            etaPhi_sqrd_prop[2] * eps_etaPhi_sqrd + 1,
            etaPhi_sqrd_prop[3] * eps_etaPhi_sqrd + 1
          )
        C2 <-
          c(etaPhi_sqrd_t[1] * eps_etaPhi_sqrd + 1,
            etaPhi_sqrd_t[2] * eps_etaPhi_sqrd + 1,
            etaPhi_sqrd_t[3] * eps_etaPhi_sqrd + 1)
        log_acceptance <-
          (
            Log_post_etaPhi(p_next, sqrt(etaPhi_sqrd_prop), sqrt(sigma_sqrd_t), mu_t, Alpha, Y) + ddirichlet(matrix(etaPhi_sqrd_t, nrow = 1), C1 , log = TRUE)
          ) - (Log_post_etaPhi(p_next, sqrt(etaPhi_sqrd_t), sqrt(sigma_sqrd_t), mu_t, Alpha, Y) + ddirichlet(matrix(etaPhi_sqrd_prop, nrow = 1), C2, log = TRUE))
        acc <- min(1, exp(log_acceptance))
        # Decision on the next step
        U = runif(1)
        if (s <- (acc>U)) etaPhi_sqrd_adapt[j+1,] <- etaPhi_sqrd_prop
        else etaPhi_sqrd_adapt[j+1,] <- etaPhi_sqrd_t
        
        ac_rate_window[j %% ac_rate_span + 1] <- s
        
        if (j >= ac_rate_span){
          current_ac_rate <- sum(ac_rate_window) / ac_rate_span
          
          # Adjust the eps_p 
          if (current_ac_rate < ac_rate_lower) eps_etaPhi_sqrd <- eps_etaPhi_sqrd*1.1
          else if (current_ac_rate > ac_rate_upper) eps_etaPhi_sqrd <- max(0, eps_etaPhi_sqrd*0.9)
        }
      }
      adapt_eps_etaPhi_sqrd <- FALSE
    }
    
    
    etaPhi_sqrd[i+1,] <- etaPhi_sqrd_next
    # Update sigma1/2
    etaPhi <- sqrt(etaPhi_sqrd[i+1,])
    sigma <- sqrt(sigma_sqrd_t)
    mu <- mu_t
    sigma1[i+1] <- sigma * (etaPhi[1] / sqrt(p_next))
    sigma2[i+1] <- sigma * (etaPhi[2] / sqrt(1 - p_next))
    mu1[i + 1] <- mu - sigma * etaPhi[3] * (sqrt(1 - p_next) / sqrt(p_next))
    mu2[i + 1] <- mu + sigma * etaPhi[3] * (sqrt(p_next) / sqrt(1 - p_next))
    
    
    ##======================Sample standard deviation========================>
    sigma_sqrd_prop <- rinvgamma(1, 2+ ((sigma_sqrd_t^2)/eps_sigma_sqrd), sigma_sqrd_t * ( 1 + (sigma_sqrd_t^2)/eps_sigma_sqrd ))
    
    
    log_acceptance <-
      (
        Log_post_sigma(p_next,sqrt(etaPhi_sqrd_next), sqrt(sigma_sqrd_prop), mu_t, ksi, Y) + dinvgamma(sigma_sqrd_t, 2+ ((sigma_sqrd_t^2)/eps_sigma_sqrd), sigma_sqrd_t * ( 1 + (sigma_sqrd_t^2)/eps_sigma_sqrd ), log = TRUE)
      ) - (Log_post_sigma(p_next, sqrt(etaPhi_sqrd_next), sqrt(sigma_sqrd_t), mu_t, ksi, Y) + dinvgamma(sigma_sqrd_prop,2+ ((sigma_sqrd_prop^2)/eps_sigma_sqrd), sigma_sqrd_prop * ( 1 + (sigma_sqrd_prop^2)/eps_sigma_sqrd ),  log = TRUE))
    acc_sigma_sqrd <- min(1, exp(log_acceptance))
    
    
    # Decision on the next step
    U = runif(1)
    if (s <- (acc_sigma_sqrd > U))
      sigma_sqrd_next <- sigma_sqrd_prop
    else
      sigma_sqrd_next <- sigma_sqrd_t
    
    # Acceptance rate
    if (i==1) ac_sigma_sqrd[i] <- s
    else ac_sigma_sqrd[i] <- ac_sigma_sqrd[i-1]+s
    ac_rate_sigma_sqrd[i] <- ac_sigma_sqrd[i]/i
    
    sigma_sqrt <- sqrt(sigma_sqrd_next)
    # Update 
    sigma1[i+1] <- sigma_sqrt * (etaPhi[1] / sqrt(p_next))
    sigma2[i+1] <- sigma_sqrt * (etaPhi[2] / sqrt(1 - p_next))
    mu1[i + 1] <- mu - sigma_sqrt * etaPhi[3] * (sqrt(1 - p_next) / sqrt(p_next))
    mu2[i + 1] <- mu + sigma_sqrt * etaPhi[3] * (sqrt(p_next) / sqrt(1 - p_next))
    
    
    if ((i%%cur_adapt==0) && ((median(ac_rate_sigma_sqrd[(i-cur_adapt/2):i])>.4) || (median(ac_rate_sigma_sqrd[(i-cur_adapt/2):i])<.3))) {
      nb_adapt <- c(nb_adapt, i)
      adapt_eps_sigma_sqrd <- TRUE
      ac_rate_window <- numeric(ac_rate_span)
    }
    
    if (adapt_eps_sigma_sqrd == TRUE){
      sigma_sqrd_adapt[1] <- sigma_sqrd_next
      for (j in 1:T_adapt){
        sigma_sqrd_t <- sigma_sqrd_adapt[j]
        sigma_sqrd_prop <- rinvgamma(1, 2+ ((sigma_sqrd_t^2)/eps_sigma_sqrd), sigma_sqrd_t * ( 1 + (sigma_sqrd_t^2)/eps_sigma_sqrd ))
        print
        # Computing acceptance probability
        log_acceptance <-
          (
            Log_post_sigma(p_next,sqrt(etaPhi_sqrd_next), sqrt(sigma_sqrd_prop), mu_t, ksi, Y) + dinvgamma(sigma_sqrd_t, 2+ ((sigma_sqrd_t^2)/eps_sigma_sqrd), sigma_sqrd_t * ( 1 + (sigma_sqrd_t^2)/eps_sigma_sqrd ), log = TRUE)
          ) - (Log_post_sigma(p_next, sqrt(etaPhi_sqrd_next), sqrt(sigma_sqrd_t), mu_t, ksi, Y) + dinvgamma(sigma_sqrd_prop,2+ ((sigma_sqrd_prop^2)/eps_sigma_sqrd), sigma_sqrd_prop * ( 1 + (sigma_sqrd_prop^2)/eps_sigma_sqrd ),  log = TRUE))
        acc_sigma_sqrd <- min(1, exp(log_acceptance))
        
        # Decision on the next step
        U = runif(1)
        if (s <- (acc_sigma_sqrd > U))
          sigma_sqrd_adapt[j+1] <- sigma_sqrd_prop
        else
          sigma_sqrd_adapt[j+1] <- sigma_sqrd_t
        
        ac_rate_window[j %% ac_rate_span + 1] <- s
        
        if (j >= ac_rate_span){
          current_ac_rate <- sum(ac_rate_window) / ac_rate_span
          
          # Adjust the eps_p
          if (current_ac_rate < ac_rate_lower) eps_sigma_sqrd <- eps_sigma_sqrd*0.95
          else if (current_ac_rate > ac_rate_upper) eps_sigma_sqrd <- eps_sigma_sqrd*1.05
        }
      }
      adapt_eps_sigma_sqrd <- FALSE
      
      
    }
    sigma <- sqrt(sigma_sqrd_next)
    etaPhi <- sqrt(etaPhi_sqrd_next)
    sigma1[i+1] <- sigma * (etaPhi[1] / sqrt(p_next))
    sigma2[i+1] <- sigma * (etaPhi[2] / sqrt(1 - p_next))
    mu1[i + 1] <- mu - sigma * etaPhi[3] * (sqrt(1 - p_next) / sqrt(p_next))
    mu2[i + 1] <- mu + sigma * etaPhi[3] * (sqrt(p_next) / sqrt(1 - p_next))
    
    sigma_sqrd[i+1] <- sigma_sqrd_next
    
    ##======================Sample mean========================>
    mu_prop <- rnorm(1, mean = mu_t, eps_mu)
    
    
    log_acceptance <-
      (
        Log_post_mu(p_next,sqrt(etaPhi_sqrd_next), sqrt(sigma_sqrd_next), mu_prop, psi, Y) + dnorm(mu_t,mu_prop, eps_mu, log = TRUE)
      ) - (Log_post_mu(p_next,sqrt(etaPhi_sqrd_next), sqrt(sigma_sqrd_next), mu_t, psi, Y) + dnorm(mu_prop,mu_t, eps_mu, log = TRUE))
    acc_mu <- min(1, exp(log_acceptance))
    
    
    # Decision on the next step
    U = runif(1)
    if (s <- (acc_mu > U))
      mu_next <- mu_prop
    else
      mu_next <- mu_t
    
    # Acceptance rate
    if (i==1) ac_mu[i] <- s
    else ac_mu[i] <- ac_mu[i-1]+s
    ac_rate_mu[i] <- ac_mu[i]/i
    
    mu <- mu_next
    # Update 
    sigma_sqrt <- sqrt(sigma_sqrd_next)
    sigma1[i+1] <- sigma_sqrt * (etaPhi[1] / sqrt(p_next))
    sigma2[i+1] <- sigma_sqrt * (etaPhi[2] / sqrt(1 - p_next))
    mu1[i+1] <- mu - sigma_sqrt * etaPhi[3] * (sqrt(1 - p_next) / sqrt(p_next))
    mu2[i+1] <- mu + sigma_sqrt * etaPhi[3] * (sqrt(p_next) / sqrt(1 - p_next))
    
    
    if ((i%%cur_adapt==0) && ((median(ac_rate_mu[(i-cur_adapt/2):i])>.4) || (median(ac_rate_mu[(i-cur_adapt/2):i])<.3))) {
      nb_adapt <- c(nb_adapt, i)
      adapt_eps_mu <- TRUE
      ac_rate_window <- numeric(ac_rate_span)
    }
    
    if (adapt_eps_mu == TRUE){
      Mu_adapt[1] <- mu_next
      for (j in 1:T_adapt){
        mu_t <- Mu_adapt[j]
        mu_prop <- rnorm(1, mean = mu_t, eps_mu)
        
        # Computing acceptance probability
        log_acceptance <-
          (
            Log_post_mu(p_next,sqrt(etaPhi_sqrd_next), sqrt(sigma_sqrd_next), mu_prop, psi, Y) + dnorm(mu_t,mu_prop, eps_mu, log = TRUE)
          ) - (Log_post_mu(p_next,sqrt(etaPhi_sqrd_next), sqrt(sigma_sqrd_next), mu_t, psi, Y) + dnorm(mu_prop,mu_t, eps_mu, log = TRUE))
        acc_mu <- min(1, exp(log_acceptance))
        # Decision on the next step
        U = runif(1)
        if (s <- (acc_mu > U))
          Mu_adapt[j+1] <- mu_prop
        else
          Mu_adapt[j+1] <- mu_t
        
        ac_rate_window[j %% ac_rate_span + 1] <- s
        
        if (j >= ac_rate_span){
          current_ac_rate <- sum(ac_rate_window) / ac_rate_span
          
          # Adjust the eps_p
          if (current_ac_rate < ac_rate_lower) eps_mu <- eps_mu*0.95
          else if (current_ac_rate > ac_rate_upper) eps_mu <- eps_mu*1.05
        }
      }
      adapt_eps_mu <- FALSE
      
      
    }
    mu <- mu_next
    sigma <- sqrt(sigma_sqrd_next)
    etaPhi <- sqrt(etaPhi_sqrd_next)
    sigma1[i+1] <- sigma * (etaPhi[1] / sqrt(p_next))
    sigma2[i+1] <- sigma * (etaPhi[2] / sqrt(1 - p_next))
    mu1[i + 1] <- mu - sigma * etaPhi[3] * (sqrt(1 - p_next) / sqrt(p_next))
    mu2[i + 1] <- mu + sigma * etaPhi[3] * (sqrt(p_next) / sqrt(1 - p_next))
    
    Mu[i+1] <- mu_next
  }
  
  n_trunc <- 100
  sigma_sqrd <- sigma_sqrd[n_trunc:length(sigma_sqrd)]
  sigma_ <- sqrt(sigma_sqrd)
  Mu <- Mu[n_trunc:length(Mu)]
  P <- P[n_trunc:length(P)]
  sigma1 <- sigma1[n_trunc:length(sigma1)]
  sigma2 <- sigma2[n_trunc:length(sigma2)]
  mu1 <- mu1[n_trunc:length(mu1)]
  mu2 <- mu2[n_trunc:length(mu2)]
  
  
  return(list(
    "mixture_wight" = P,
    "eta1" = sqrt(etaPhi_sqrd[,1]),
    "eta2" = sqrt(etaPhi_sqrd[,2]),
    "phi"  = sqrt(etaPhi_sqrd[,3]),
    "mu" = Mu,
    "sigma" = sigma_,
    "mu1" = mu1,
    "mu2" = mu2,
    "sigma1" = sigma1,
    "sigma2" = sigma2
  ))
}
