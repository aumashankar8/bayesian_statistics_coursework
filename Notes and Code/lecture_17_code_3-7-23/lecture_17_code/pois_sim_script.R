###
###  Simulate Count Regression Data
###
setwd("~/UT Austin/Spring 2023/Bayesian Stats/Notes and Code/lecture_17_code/lecture_17_code")
n=30
p=3
X=matrix(1,n,p)
set.seed(101)
X[,2]=rnorm(n)
X[,3]=rnorm(n)
beta=c(.5,-1,2)
lambda=exp(X%*%beta)
y=rpois(n,lambda)  # simulate data

pairs(cbind(y,X[,-1]))
summary(glm(y~0+X,family=poisson(link="log")))  # check MLE 

###
###  Specify Count Regression Model in JAGS (FINISH THIS PART)
###

library(rjags)

m.jags <-"
  model{
    for(i in 1:n){
      y[i] ~ dpois(lambda[i])
      log(lambda[i]) = b0 + b1*X[i,2] + b2*X[i,3]
    }
    b0 ~ dnorm(mu,tau)
    b1 ~ dnorm(mu,tau)
    b2 ~ dnorm(mu,tau)
}
"

mod<-textConnection(m.jags)
###
###  Fit Count Regression Model w/ JAGS 
###

mu=0
tau=1/100
m.out<-jags.model(mod,data=list('y'=y,'n'=n,'X'=X,'mu'=mu,'tau'=tau),n.chains=1,n.adapt=0) # build model and algorithm

n.mcmc=10000
n.burn=round(.2*n.mcmc)
update(m.out,n.burn) # perform burn-in
m.samples=jags.samples(m.out,c('b0','b1','b2'),n.mcmc) # fit model post-burn-in

beta.post.mat=cbind(m.samples$b0[1,,1],m.samples$b1[1,,1],m.samples$b2[1,,1])
matplot(beta.post.mat,type="l",lty=1,ylab=bquote(beta))
abline(h=beta,col=8)

###
###  Count Regression Inference 
###

library(vioplot)
vioplot(data.frame(beta.post.mat),names=expression(beta[0],beta[1],beta[2]))
abline(h=0,col=8)

apply(beta.post.mat,2,mean) # marginal posterior means for beta 
apply(beta.post.mat,2,sd) # marginal posterior means for beta 
apply(beta.post.mat,2,quantile,c(0.025,.975)) # marginal posterior 95% CI for beta 


