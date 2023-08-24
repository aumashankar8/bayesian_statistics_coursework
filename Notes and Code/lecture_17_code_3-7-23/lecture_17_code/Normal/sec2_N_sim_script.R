###
###  Simulate Linear Regression Data
###

n=30
p=3
X=matrix(1,n,p)
set.seed(101)
X[,2]=runif(n)
X[,3]=runif(n)
beta=c(1,-.5,2)
sig=.1
y=rnorm(n,X%*%beta,sig)  # simulate data

pairs(cbind(y,X[,-1]))
summary(lm(y~0+X))  # check least squares fit 

###
###  Specify Linear Regression Model in JAGS
###

library(rjags)

m.jags <-"
  model{
    for(i in 1:n){
      y[i] ~ dnorm(mu[i], tau)  # note: tau is precision=1/variance!  
      mu[i] = X[i,1:p]%*%beta
    }
    beta ~ dmnorm(mu.beta,Tau.beta)
    tau <- pow(sig, -2)
    sig ~ dunif(0,10)
}
"

mod<-textConnection(m.jags)

###
###  Fit Linear Regression Model w/ JAGS 
###

mu.beta=rep(0,p)
Sig.beta=100*diag(p)
Tau.beta=solve(Sig.beta)
m.out<-jags.model(mod,data=list('y'=y,'n'=n,'X'=X,'p'=p,'mu.beta'=mu.beta,'Tau.beta'=Tau.beta),n.chains=1,n.adapt=0) # build model and algorithm

n.mcmc=10000
n.burn=round(.2*n.mcmc)
update(m.out,n.burn) # perform burn-in
m.samples=jags.samples(m.out,c('beta','sig'),n.mcmc) # fit model post-burn-in

layout(matrix(1:2,2,1))
matplot(t(m.samples$beta[,,1]),type="l",lty=1,ylab=bquote(beta))
abline(h=beta,col=8)
plot(m.samples$sig[1,,1],type="l",lty=1,ylab=bquote(sigma))
abline(h=sig,col=8)

###
###  Linear Regression Inference 
###

library(vioplot)
vioplot(data.frame(t(m.samples$beta[,,1])),names=expression(beta[0],beta[1],beta[2]))
abline(h=0,col=8)

apply(m.samples$beta[,,1],1,mean) # marginal posterior means for beta 
apply(m.samples$beta[,,1],1,sd) # marginal posterior means for beta 
apply(m.samples$beta[,,1],1,quantile,c(0.025,.975)) # marginal posterior 95% CI for beta 

mean(m.samples$sig[1,,1])  # marginal posterior mean for sigma
sd(m.samples$sig[1,,1])  # marginal posterior mean for sigma

###
###  Watch Out for Transformations!
###

mean(m.samples$sig[1,,1]^2)  # marginal posterior mean for sigma^2
(mean(m.samples$sig[1,,1])^2)  # This is incorrect!  



