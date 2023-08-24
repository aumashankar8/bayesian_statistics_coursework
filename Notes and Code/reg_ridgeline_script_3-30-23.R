###
###  Simulate Linear Regression Data
###

set.seed(30003)

n=50
p=10
X=matrix(1,n,p)
set.seed(101)
for(j in 2:p){
  X[,j]=runif(n)
}
beta=runif(p,-2,2)
sig=.1
y=rnorm(n,X%*%beta,sig)  # simulate data

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
m.out<-jags.model(mod,data=list('y'=y,'n'=n,'X'=X,'p'=p,'mu.beta'=mu.beta,'Tau.beta'=Tau.beta),n.chains=3,n.adapt=0) # build model and algorithm

n.mcmc=10000
n.burn=round(.2*n.mcmc)
update(m.out,n.burn) # perform burn-in
m.samples=coda.samples(m.out,c('beta','sig'),n.mcmc) # fit model post-burn-in

gelman.diag(m.samples)

layout(matrix(1:2,2,1))
matplot(m.samples[[1]][,1:p],type="l",lty=1,ylab=bquote(beta))
abline(h=beta,col=8)
plot(c(m.samples[[1]][,'sig']),type="l",lty=1,ylab=bquote(sigma))
abline(h=sig,col=8)

###
###  Make Ridgeline Plot for Slope Parameters 
###

library(forcats)
library(ggridges)
library(ggplot2)

dist.labels=as.factor(rep(as.character(1:(p-1)),each=n.mcmc))
slopes.mat=as.matrix(m.samples[[1]][,2:p])
slopes.df=data.frame(x=c(slopes.mat),y=dist.labels)
ggplot(slopes.df, aes(x=x,y=fct_reorder(dist.labels,rep(2:p,each=n.mcmc))))+geom_density_ridges(rel_min_height = 0.005)+labs(y="posterior density",x="parameter value",title="Slope Coefficients")+geom_vline(xintercept=0,size=1,linetype="dashed",col=rgb(0,0,0,.5))

