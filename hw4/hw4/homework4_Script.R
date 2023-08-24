###   Homework 4

###   Libraries
library(rjags)
library(vioplot)
###
###   Read the Data
###


 
setwd("~/UT Austin/Spring 2023/Bayesian Stats/hw4/hw4")
mosquito=read.table("mosquitofish.txt",header=TRUE) # read in mosquitofish data
head(mosquito)


n  = 32 # those groups of fish that were only counted twice
y  = mosquito$N1[1:n]
N  = mosquito$N0[1:n]
l1 = mosquito$L1[1:n]
d  = mosquito$d[1:n]

l1_stand = scale(l1)[,1] #standardized L1 - mean length of surviving fish
d_stand = scale(d)[,1] #standardized Days of study

l1_d_product = l1_stand * d_stand


## Jags Regression
# (a)

# Covariates - Intercept, stand l1, stand d, product of l1 and d
p = 4
X = matrix(1,n,p)
X[,2] = l1_stand
X[,3] = d_stand
X[,4] = l1_d_product

mosquito_jags <-"
  model{
    for(i in 1:n){
      y[i] ~ dbin(theta[i],N[i]) 
      probit(theta[i]) = b0 + b1*X[i,2] + b2*X[i,3] + b3*(X[i,2]*X[i,3])
    }
    b0 ~ dnorm(mu,tau)
    b1 ~ dnorm(mu,tau)
    b2 ~ dnorm(mu,tau)
    b3 ~ dnorm(mu,tau)
}
"
mod <- textConnection(mosquito_jags)
# 
# mu.beta=rep(0,p)
# Sig.beta=100*diag(p)
# Tau.beta=solve(Sig.beta)

mu=0
tau=1/2.25 #Precision. 1 / sigma^2
m.out<-jags.model(mod,data=list('y'=y,'n'=n, 'N'=N, 'X'=X,'mu'=mu,'tau'=tau),n.chains=1,n.adapt=0) # build model and algorithm

n.mcmc=10000
n.burn=round(.2*n.mcmc)
update(m.out,n.burn) # perform burn-in
m.samples=jags.samples(m.out,c('b0','b1','b2','b3'),n.mcmc) # fit model post-burn-in

#Trace Plot for (a)
beta.post.mat=cbind(m.samples$b0[1,,1],m.samples$b1[1,,1],m.samples$b2[1,,1],m.samples$b3[1,,1])
matplot(beta.post.mat,type="l",lty=1,ylab=bquote(beta))

#Violin plot for (a)
vioplot(data.frame(beta.post.mat),names=expression(beta[0],beta[1],beta[2],beta[3]))
abline(h=0,col=8)

apply(beta.post.mat,2,mean) # marginal posterior means for beta 
apply(beta.post.mat,2,sd) # marginal posterior sd for beta 
apply(beta.post.mat,2,quantile,c(0.025,.975)) # marginal posterior 95% CI for beta 


# Covariates - Intercept, stand l1, stand d
# (b)

mosquito_jags_b <-"
  model{
    for(i in 1:n){
      y[i] ~ dbin(theta[i],N[i]) 
      probit(theta[i]) = b0 + b1*X[i,2] + b2*X[i,3]
    }
    b0 ~ dnorm(mu,tau)
    b1 ~ dnorm(mu,tau)
    b2 ~ dnorm(mu,tau)
}
"
mod_b <- textConnection(mosquito_jags_b)
# 
# mu.beta=rep(0,p)
# Sig.beta=100*diag(p)
# Tau.beta=solve(Sig.beta)

mu=0
tau=1/2.25 #Precision. 1 / sigma^2
m.out_b<-jags.model(mod_b,data=list('y'=y,'n'=n, 'N'=N, 'X'=X,'mu'=mu,'tau'=tau),n.chains=1,n.adapt=0) # build model and algorithm

n.mcmc=10000
n.burn=round(.2*n.mcmc)
update(m.out_b,n.burn) # perform burn-in
m.samples_b=jags.samples(m.out_b,c('b0','b1','b2'),n.mcmc) # fit model post-burn-in

#Trace Plot for (b)
beta.post.mat_b=cbind(m.samples_b$b0[1,,1],m.samples_b$b1[1,,1],m.samples_b$b2[1,,1])
matplot(beta.post.mat_b,type="l",lty=1,ylab=bquote(beta))

#Violin plot for (b)
vioplot(data.frame(beta.post.mat_b),names=expression(beta[0],beta[1],beta[2])) 
abline(h=0,col=8)

apply(beta.post.mat_b,2,mean) # marginal posterior means for beta 
apply(beta.post.mat_b,2,sd) # marginal posterior sd for beta 
apply(beta.post.mat_b,2,quantile,c(0.025,.975)) # marginal posterior 95% CI for beta 

# Covariates - Intercept, stand l1
# (c)

mosquito_jags_c <-"
  model{
    for(i in 1:n){
      y[i] ~ dbin(theta[i],N[i]) 
      probit(theta[i]) = b0 + b1*X[i,2]
    }
    b0 ~ dnorm(mu,tau)
    b1 ~ dnorm(mu,tau)
}
"
mod_c <- textConnection(mosquito_jags_c)
# 
# mu.beta=rep(0,p)
# Sig.beta=100*diag(p)
# Tau.beta=solve(Sig.beta)

mu=0
tau=1/2.25 #Precision. 1 / sigma^2
m.out_c<-jags.model(mod_c,data=list('y'=y,'n'=n, 'N'=N, 'X'=X,'mu'=mu,'tau'=tau),n.chains=1,n.adapt=0) # build model and algorithm

n.mcmc=10000
n.burn=round(.2*n.mcmc)
update(m.out_c,n.burn) # perform burn-in
m.samples_c=jags.samples(m.out_c,c('b0','b1'),n.mcmc) # fit model post-burn-in

#Trace Plot for (c)
beta.post.mat_c=cbind(m.samples_c$b0[1,,1],m.samples_c$b1[1,,1],m.samples_c$b2[1,,1])
matplot(beta.post.mat_c,type="l",lty=1,ylab=bquote(beta))

#Violin plot for (c)
vioplot(data.frame(beta.post.mat_c),names=expression(beta[0],beta[1])) 
abline(h=0,col=8)

apply(beta.post.mat_c,2,mean) # marginal posterior means for beta 
apply(beta.post.mat_c,2,sd) # marginal posterior sd for beta 
apply(beta.post.mat_c,2,quantile,c(0.025,.975)) # marginal posterior 95% CI for beta 


# Covariates - Intercept, stand d
# (d)

mosquito_jags_d <-"
  model{
    for(i in 1:n){
      y[i] ~ dbin(theta[i],N[i]) 
      probit(theta[i]) = b0 + b1*X[i,3]
    }
    b0 ~ dnorm(mu,tau)
    b1 ~ dnorm(mu,tau)
}
"
mod_d <- textConnection(mosquito_jags_d)
# 
# mu.beta=rep(0,p)
# Sig.beta=100*diag(p)
# Tau.beta=solve(Sig.beta)

mu=0
tau=1/2.25 #Precision. 1 / sigma^2
m.out_d<-jags.model(mod_d,data=list('y'=y,'n'=n, 'N'=N, 'X'=X,'mu'=mu,'tau'=tau),n.chains=1,n.adapt=0) # build model and algorithm

n.mcmc=10000
n.burn=round(.2*n.mcmc)
update(m.out_d,n.burn) # perform burn-in
m.samples_d=jags.samples(m.out_d,c('b0','b1'),n.mcmc) # fit model post-burn-in

#Trace Plot for (d)
beta.post.mat_d=cbind(m.samples_d$b0[1,,1],m.samples_d$b1[1,,1])
matplot(beta.post.mat_d,type="l",lty=1,ylab=bquote(beta))

#Violin plot for (d)
vioplot(data.frame(beta.post.mat_d),names=expression(beta[0],beta[1])) 
abline(h=0,col=8)

apply(beta.post.mat_d,2,mean) # marginal posterior means for beta 
apply(beta.post.mat_d,2,sd) # marginal posterior sd for beta 
apply(beta.post.mat_d,2,quantile,c(0.025,.975)) # marginal posterior 95% CI for beta 

# Covariates - Intercept ONLY
# (e)

mosquito_jags_e <-"
  model{
    for(i in 1:n){
      y[i] ~ dbin(theta[i],N[i]) 
      probit(theta[i]) = b0
    }
    b0 ~ dnorm(mu,tau)
}
"
mod_e <- textConnection(mosquito_jags_e)
# 
# mu.beta=rep(0,p)
# Sig.beta=100*diag(p)
# Tau.beta=solve(Sig.beta)

mu=0
tau=1/2.25 #Precision. 1 / sigma^2
m.out_e<-jags.model(mod_e,data=list('y'=y,'n'=n, 'N'=N, 'X'=X,'mu'=mu,'tau'=tau),n.chains=1,n.adapt=0) # build model and algorithm

n.mcmc=10000
n.burn=round(.2*n.mcmc)
update(m.out_e,n.burn) # perform burn-in
m.samples_e=jags.samples(m.out_e,c('b0'),n.mcmc) # fit model post-burn-in

#Trace Plot for (e)
beta.post.mat_e=cbind(m.samples_e$b0[1,,1])
matplot(beta.post.mat_e,type="l",lty=1,ylab=bquote(beta))

#Violin plot for (e)
vioplot(data.frame(beta.post.mat_e),names=expression(beta[0])) 
abline(h=0,col=8)

apply(beta.post.mat_e,2,mean) # marginal posterior means for beta 
apply(beta.post.mat_e,2,sd) # marginal posterior sd for beta 
apply(beta.post.mat_e,2,quantile,c(0.025,.975)) # marginal posterior 95% CI for beta 

#-----------------------------------------------------------------------------------#
### PLOTS and Quantiles

## All Trace Plots 
matplot(beta.post.mat,type="l",lty=1,ylab=bquote(beta))   #(a)
matplot(beta.post.mat_b,type="l",lty=1,ylab=bquote(beta)) #(b)
matplot(beta.post.mat_c,type="l",lty=1,ylab=bquote(beta)) #(c)
matplot(beta.post.mat_d,type="l",lty=1,ylab=bquote(beta)) #(d)
matplot(beta.post.mat_e,type="l",lty=1,ylab=bquote(beta)) #(e)

## All Violin Plots

vioplot(data.frame(beta.post.mat),names=expression(beta[0],beta[1],beta[2],beta[3]))  #(a)
abline(h=0,col=8)

vioplot(data.frame(beta.post.mat_b),names=expression(beta[0],beta[1],beta[2]))        #(b)
abline(h=0,col=8)

vioplot(data.frame(beta.post.mat_c),names=expression(beta[0],beta[1]))                #(c)
abline(h=0,col=8)

vioplot(data.frame(beta.post.mat_d),names=expression(beta[0],beta[1]))                #(d)
abline(h=0,col=8)

vioplot(data.frame(beta.post.mat_e),names=expression(beta[0]))                        #(e)
abline(h=0,col=8)

## All Means
apply(beta.post.mat,2,mean)   #(a)
apply(beta.post.mat_b,2,mean) #(b)
apply(beta.post.mat_c,2,mean) #(c)
apply(beta.post.mat_d,2,mean) #(d)
apply(beta.post.mat_e,2,mean) #(e)


## All Quantiles
apply(beta.post.mat,2,quantile,c(0.025,.975))   #(a)
apply(beta.post.mat_b,2,quantile,c(0.025,.975)) #(b)
apply(beta.post.mat_c,2,quantile,c(0.025,.975)) #(c)
apply(beta.post.mat_d,2,quantile,c(0.025,.975)) #(d)
apply(beta.post.mat_e,2,quantile,c(0.025,.975)) #(e)


## DIC Calculations

# (a) DIC = 74.6349860620726
# (b) DIC = 97.1322909927127
# (c) DIC = 136.613174557372
# (d) DIC = 137.52021356181
# (e) DIC = 139.469354229737

#(a)

postbeta = apply(beta.post.mat, 2, mean)
posttheta = pnorm(X%*%postbeta)
Dhat = mean(-2*sum(log(dbinom(y,N,posttheta))))
pD = 2 * sum(log(mean(exp(beta.post.mat), 2)))
DIC=Dhat+2*pD
print(paste("DIC =",DIC))
# "DIC = 74.6349860620726"

#(b)
postbeta = apply(beta.post.mat_b, 2, mean)
posttheta = pnorm(X[,-4]%*%postbeta)
Dhat = mean(-2*sum(log(dbinom(y,N,posttheta))))
pD = 2 * sum(log(mean(exp(beta.post.mat), 2)))
DIC=Dhat+2*pD
print(paste("DIC =",DIC))
# "DIC = 97.1322909927127"

#(c)
postbeta = apply(beta.post.mat_c, 2, mean)
posttheta = pnorm(X[,1:2]%*%postbeta)
Dhat = mean(-2*sum(log(dbinom(y,N,posttheta))))
pD = 2 * sum(log(mean(exp(beta.post.mat), 2)))
DIC=Dhat+2*pD
print(paste("DIC =",DIC))
# "DIC = 136.613174557372"

#(d)
postbeta = apply(beta.post.mat_d, 2, mean)
posttheta = pnorm(X[,c(1,3)]%*%postbeta)
Dhat = mean(-2*sum(log(dbinom(y,N,posttheta))))
pD = 2 * sum(log(mean(exp(beta.post.mat), 2)))
DIC=Dhat+2*pD
print(paste("DIC =",DIC))
# "DIC = 137.52021356181"

#(e)
postbeta = apply(beta.post.mat_e, 2, mean)
posttheta = pnorm(X[,1]*postbeta)
Dhat = mean(-2*sum(log(dbinom(y,N,posttheta))))
pD = 2 * sum(log(mean(exp(beta.post.mat), 2)))
DIC=Dhat+2*pD
print(paste("DIC =",DIC))
# "DIC = 139.469354229737"

