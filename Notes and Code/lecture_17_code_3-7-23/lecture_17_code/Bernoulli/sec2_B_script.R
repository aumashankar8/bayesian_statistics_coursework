###
###  Read in Coexistence Data
###
#
#  Hein, C.L., Öhlund, G. and Englund, G., 2014. Fish introductions reveal the temperature dependence of species interactions. Proceedings of the Royal Society B: Biological Sciences, 281(1775), p.20132641.
#
setwd("~/UT Austin/Spring 2023/Bayesian Stats/Notes and Code/lecture_17_code/lecture_17_code/Bernoulli")
logit<-function(theta){log(theta)-log(1-theta)}
logit.inv<-function(x){exp(x)/(1+exp(x))}

ce.df=read.table("coexist.txt",header=TRUE)
head(ce.df)

n=dim(ce.df)[1]
y=ce.df$coexist
p=4

X=matrix(1,n,p)
X[,2]=scale(ce.df$elev)
X[,3]=scale(ce.df$maxdepth)
X[,4]=scale(ce.df$temp1)

pairs(cbind(y,X[,-1]))
summary(glm(y~0+X,family="binomial"))  # check MLE 

###
###  Specify Binary Regression Model in JAGS
###

library(rjags)

m.jags <-"
  model{
    for(i in 1:n){
      y[i] ~ dbin(theta[i],N) 
      logit(theta[i]) <- b0+b1*X[i,2]+b2*X[i,3]+b3*X[i,4]
    }
    b0 ~ dnorm(mu,tau)
    b1 ~ dnorm(mu,tau)
    b2 ~ dnorm(mu,tau)
    b3 ~ dnorm(mu,tau)
}
"

mod<-textConnection(m.jags)

###
###  Fit Binary Regression Model w/ JAGS 
###

mu=0
tau=1/2.25 #Precision. 1 / sigma^2
m.out<-jags.model(mod,data=list('y'=y,'n'=n,'X'=X,'N'=1,'mu'=mu,'tau'=tau),n.chains=1,n.adapt=0) # build model and algorithm

n.mcmc=50000
n.burn=round(.2*n.mcmc)
update(m.out,n.burn) # perform burn-in
m.samples=jags.samples(m.out,c('b0','b1','b2','b3'),n.mcmc) # fit model post-burn-in

beta.post.mat=cbind(m.samples$b0[1,,1],m.samples$b1[1,,1],m.samples$b2[1,,1],m.samples$b3[1,,1])
matplot(beta.post.mat,type="l",lty=1,ylab=bquote(beta))

###
###  Binary Regression Inference 
###

library(vioplot) #install this)
vioplot(data.frame(beta.post.mat),names=expression(beta[0],beta[1],beta[2],beta[3]))
abline(h=0,col=8)

apply(beta.post.mat,2,mean) # marginal posterior means for beta 
apply(beta.post.mat,2,sd) # marginal posterior sd for beta 
apply(beta.post.mat,2,quantile,c(0.025,.975)) # marginal posterior 95% CI for beta 


