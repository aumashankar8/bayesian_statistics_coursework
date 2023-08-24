###
###  Read in Data for Binomial Model
###
setwd("~/UT Austin/Spring 2023/Bayesian Stats/Notes and Code/lecture_16_code_3-2-23/lecture_16_code")
mf.df=read.table("mosquitofish.txt",header=TRUE) # read in mosquitofish data
head(mf.df)

n=32 # those groups of fish that were only counted twice
y=mf.df$N1[1:n]
N=mf.df$N0[1:n]

plot(N,y,asp=TRUE)
abline(0,1)

###
###  Fit Binomial Model using 1 Chain
###

library(rjags)
#dbin NOT dbinom. dbin also has different inputs so be careful
#JAGS uses stuff from antiquated syntax
m.jags <-"
  model{
    for(i in 1:n){
      y[i] ~ dbin(p,N[i]) 
    }
    p ~ dbeta(alpha,beta)
}
"

n.mcmc=1000
n.burn=round(.2*n.mcmc)
mod<-textConnection(m.jags) # read model, do the model statement in another script 
m.out<-jags.model(mod,data=list('y'=y,'n'=n,'N'=N,'alpha'=1,'beta'=1),n.chains=1,n.adapt=0) # build model graph | n.cahins is amount of models it will produce and then converge | n.adapt is number of adaptations to adapt. n
m.samples=jags.samples(m.out,c('p'),n.mcmc) # fit model to data
#jags.samples(model, which var to write out, how many mcmc iterations)

#Output. p parameter[which p,which mcmc iters (all of em), which chain]
plot(m.samples$p[1,,1],type="l",lty=1,ylab="p",xlab="iteration",ylim=c(0,1)) # trace plot

hist(m.samples$p[1,n.burn:n.mcmc,1],breaks=30,xlim=c(0,1),prob=TRUE,xlab="p",ylab="density",main="") # posterior histogram 
mean(m.samples$p[1,n.burn:n.mcmc,1]) #Posterior mean = 0.874

lines(density(m.samples$p[1,n.burn:n.mcmc,1]),lwd=2,col=2) 

###
###  Fit Binomial Model using 3 Chains and calculate R-hat by hand (Gelman et al., 2014; BDA3)
###

library(rjags)

m.jags <-"
  model{
    for(i in 1:n){
      y[i] ~ dbin(p,N[i])
    }
    p ~ dbeta(alpha,beta)
}
"

n.mcmc=10000
n.burn=round(.2*n.mcmc)
mod<-textConnection(m.jags) # read model
m.out<-jags.model(mod,data=list('y'=y,'n'=n,'N'=N,'alpha'=1,'beta'=1),n.chains=3,n.adapt=0) # build model graph 
m.samples=jags.samples(m.out,c('p'),n.mcmc) # fit model to data

#m.samples$p[1,,1]=rbeta(n.mcmc,8,2)

matplot(m.samples$p[1,,],type="l",col=1:3,lty=1,ylab="p",xlab="iteration",ylim=c(0,1)) # trace plots

keep.idx=n.burn:n.mcmc # discard burn-in
K=length(keep.idx)
w=mean(apply(m.samples$p[1,keep.idx,],2,var))
b=K*var(apply(m.samples$p[1,keep.idx,],2,mean))
v.hat=w*(K-1)/K+b/K
r.hat=sqrt(v.hat/w)
r.hat #When it's not great, this will be large. Want it close to 1

###
###  Fit Binomial Model using 3 Chains and use CODA to calculate R-hat (using different formulation) 
###

library(rjags)
library(coda)

m.jags <-"
  model{
    for(i in 1:n){
      y[i] ~ dbin(p,N[i])
    }
    p ~ dbeta(alpha,beta)
}
"

n.mcmc=100
n.burn=round(.2*n.mcmc)
mod<-textConnection(m.jags) # read model
m.out<-jags.model(mod,data=list('y'=y,'n'=n,'N'=N,'alpha'=1,'beta'=1),n.chains=3,n.adapt=0) # build model graph 
m.samples=coda.samples(m.out,c('p'),n.mcmc) # fit model to data
summary(m.samples) #Gives a mean
#codamenu() lets you get some extra things
    str(m.samples)
    plot(m.samples[[1]], type = "l") #Double bracket is for each chain


gelman.diag(m.samples) #point estimate is r_hat. Below 1.1 is small. 

###
###  Simulate Data for Hierarchical Binomial Model
###

n=10
lambda=20
N=rpois(n,lambda)
p=.8
y=rbinom(n,N,p)

plot(N,y,xlim=c(0,max(N)),ylim=c(0,max(N)))
abline(0,1)

###
###  Fit Hierarchical Binomial-Poisson Model assuming p and all N Unknown (but lambda known)
###
m.jags <-"
  model{
    for(i in 1:n){
      y[i] ~ dbin(p,N[i])
      N[i] ~ dpois(lambda)
    }
    p ~ dbeta(alpha,beta)
}
"
#Can try logit(p) ~ dnorm(mu, tau) , tau is inverse of sigma_sqrd

n.mcmc=10000
n.burn=round(.2*n.mcmc)
mod<-textConnection(m.jags) # read model
m.out<-jags.model(mod,data=list('y'=y,'n'=n,'alpha'=1,'beta'=1,'lambda'=lambda),n.chains=3,n.adapt=0) # build model graph 
#Jags Sample
m.samples <- jags.samples(m.out, c('p','N'), n.mcmc)
#Coda Sample
m.samples.coda <- coda.samples(m.out, c('p','N'), n.mcmc)
summary(m.samples.coda)
gelman.diag(m.samples.coda)

#Traceplots using jags.samples
matplot(m.samples$p[1,,],type="l",col=1:3,lty=1,ylab="p",xlab="iteration",ylim=c(0,1)) # trace plot for p
matplot(m.samples$N[1,,],type="l",col=1:3,lty=1,ylab="N",xlab="iteration") # trace plot for N# trace plot for N
mean(m.samples$p[1,n.burn:n.mcmc,1]) #Posterior mean = 0.874

hist(m.samples$p[1,n.burn:n.mcmc,1],breaks=30,xlim=c(0,1),prob=TRUE,xlab="p",ylab="density",main="") # posterior histogram of p
hist(m.samples$N[1,n.burn:n.mcmc,1],breaks=30,prob=TRUE,xlab="N",ylab="density",main="") # posterior histogram of N

#R_Hat for p
keep.idx=n.burn:n.mcmc # discard burn-in
K=length(keep.idx)
w=mean(apply(m.samples$p[1,keep.idx,],2,var))
b=K*var(apply(m.samples$p[1,keep.idx,],2,mean))
v.hat=w*(K-1)/K+b/K
r.hat=sqrt(v.hat/w)
r.hat 

#R_hat for N
keep.idx=n.burn:n.mcmc # discard burn-in
K=length(keep.idx)
w=mean(apply(m.samples$N[1,keep.idx,],2,var))
b=K*var(apply(m.samples$N[1,keep.idx,],2,mean))
v.hat=w*(K-1)/K+b/K
r.hat=sqrt(v.hat/w)
r.hat 
