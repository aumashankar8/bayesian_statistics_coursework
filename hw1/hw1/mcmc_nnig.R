mcmc.nnig <- function(y, mu0, sigma20, v, sigma2.start, n.mcmc){
  
  # number of data points
  n <- length(y) 
  # vector that holds mu updates
  mu.save <- rep(0, n.mcmc)
  # vector that holds sigma2 updates
  sigma2.save <- rep(0, n.mcmc)
  
  # Gibbs sampler
  
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
    
    # parameters for the full conditional distribution of sigma2
    q = (v + n)/2
    r_den <- 0.5 + sum((y-mu.k)^2)
    r <- 1/r_den
    
    # sample sigma2 from the full conditional distribution of sigma2
    sigma2.k <- 1/rgamma(1, shape=q ,scale = r)
    
    #save sigma update
    sigma2.save[k] <- sigma2.k
  } 
  
  list(mu.save=mu.save, sigma2.save=sigma2.save, n.mcmc=n.mcmc)
}