mcmc.randomWalk <- function(y, mu0, sigma20, mu_log, sig_log, sigma2.start, sig.tune, n.mcmc){
  
  # number of data points
  n <- length(y) 
  # vector that holds mu updates
  mu.save <- rep(0, n.mcmc)

  
  # Gibbs sampler
  #Implements M-H sampler for Normal-Normal, where y~N(mu,sigma2), log(sigma)~N(mu_log,sig_log)
  # initialize mu
  sigma2.k <- sigma2.start
  # Gibbs updates for k=1 to the upper limit of MCMC iterations
  for (k in 1:n.mcmc){
    # parameters for the full conditional distribution of mu
    a_num <- sigma20*sum(y) + mu0*sigma2.k
    a_den <- n*sigma20 + sigma2.k
    a <- a_num/a_den
    
    b_num <- sigma2.k*sigma20
    b_den <- n*sigma20 + sigma2.k
    b <- b_num/b_den
    
    # sample mu from full conditional posterior of mu
    mu.k <- rnorm(1, mean=a, sd=sqrt(b))
    
    # save mu update
    mu.save[k] = mu.k
  }
  #-------------------------------------------
  
  #Random Walk Proposal for log(sigma)
  library(mvtnorm)
  
  #Setup Variables
  # vector that holds sigma2 updates
  sigma2.save <- rep(0, n.mcmc)
  # vector that holds log_sigma updates
  log_sigma.save <- rep(0, n.mcmc)
  
  mh.sigma <- 1
  sigma2 <- 0
  log_sigma <- 0.5 * log(sigma20) #This transforms the sigma2 hyper parameter to log(sigma)
  ####
  ####  Begin Random Walk
  ####
  
  for(k in 1:n.mcmc){
    if(k%%1000==0)  cat(k," ")
    
    ####
    ####  Sample log(sigma) 
    ####
    
    log_sigma.star <- rnorm(1,log_sigma, sig.tune)
    sigma2.star <- exp(log_sigma.star)^2
    
    mh1 <- sum(dnorm(y, mu0, sigma2.star, log = TRUE)) + dnorm(log_sigma.star,mu_log,sig_log,log = TRUE)
    mh2 <- sum(dnorm(y, mu0, sigma20, log = TRUE)) + dnorm(log_sigma, mu_log, sig_log, log = TRUE)
    mh=exp(mh1-mh2)
    
    if(mh > runif(1)){
      sigma2 <- sigma2.star
      log_sigma <- log_sigma.star
      mh.sigma <- mh.sigma + 1
    }
    
    ####
    ####  Save Samples 
    ####
    sigma2.k <- sigma2
    log_sigma.k <- log_sigma
    
    #save sigma update
    sigma2.save[k] <- sigma2.k
    log_sigma.save[k] <- log_sigma.k
  }
  
  
  list(mu.save = mu.save, 
       sigma2.save = sigma2.save, 
       log_sigma.save = log_sigma.save,
       n.mcmc = n.mcmc)
} 