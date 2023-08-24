# set working directory
setwd("~/UT Austin/Spring 2023/Bayesian Stats/hw1/hw1")

# source MCMC for Gibbs sampling
source("mcmc_nnig.R")

# source invchisquare function
source("invchisquare.R")

# read medical data
df.med_data <- read.csv("MedicalData.csv")

# subset young adults data
df.young_adult <- df.med_data[df.med_data$Age == "young adult",]

# get y vector (BMI of young adults)
y = df.young_adult$BMI

# Hyper parameters for mean
mu0 = mean(y)
sigma20 = var(y)

# Hyper parameter for v
v = 1/var(y) + 2

# MCMC (Gibbs) for Normal-Normal-InverseChiSquare
mcmc.out <- mcmc.nnig(y, mu0, sigma20, v, sigma2.start=25, n.mcmc=1000000)

# trace plot
plot(mcmc.out$mu.save, type='l', xlab='k', ylab=bquote(mu), 
     main =bquote("Trace Plot for " ~ mu ))

plot(mcmc.out$sigma2.save, type='l', xlab='k', ylab=bquote(sigma^2),
     main =bquote("Trace Plot for " ~ sigma^2 ))

# prior distributions
# density of the prior distribution of mu
curve(dnorm(x,mu0,sqrt(sigma20)),from=0, to=60, type='l',
      xlab=bquote(mu), ylab='Density', main =bquote("[" ~ mu ~ "]"))
curve(dinvchisquare(x,v), from=0, to=2.5, type='l', xlab=bquote(sigma^2),
      ylab='Density', main =bquote("[" ~ sigma^2 ~ "]"))

# posterior distributions
# marginal posterior distribution of mu
hist(mcmc.out$mu.save, col=8, probability=TRUE, xlab = bquote(mu),
     main =bquote("[" ~ mu ~ "| \U00B7 ]")) 

# density of the prior distribution of mu
curve(dnorm(x,mu0,sqrt(sigma20)), from=15, to=50, col=8, lty=2, add=TRUE)
# marginal posterior distribution of sigma2
hist(mcmc.out$sigma2.save,xlim=c(0,300), col=8, probability=TRUE,
     xlab = bquote(sigma^2),
     main =bquote("[" ~ sigma^2 ~ "| \U00B7 ]"))
# density of the prior distribution of sigma2
curve(dinvchisquare(x,v), from=0, to=300, col=8, lty=2, add=TRUE)

# autocorrelation plots
acf(mcmc.out$mu.save, main=bquote("Autocorrelation for " ~ mu))
acf(mcmc.out$sigma2.save, main=bquote("Autocorrelation for " ~ sigma^2))

# 95% CI
CI_mean = quantile(mcmc.out$mu.save, probs = c(0.025, 0.975))
CI_var = quantile(mcmc.out$sigma2.save, probs = c(0.025, 0.975))

# calculate P(mu>30|y)
mu.posterior = mcmc.out$mu.save
prob = length(mu.posterior[mu.posterior>30])/length(mu.posterior)
print(paste('Probability of mean greater than 30 = ', prob))

