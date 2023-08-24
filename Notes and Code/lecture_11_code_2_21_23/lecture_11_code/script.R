####
####  Read in Data
####

bird.df=read.csv("birds.csv",header=TRUE)
bird.df
N=dim(bird.df)[1]

####
####  Scale Covariates 
####

bird.df$area.sc=scale(bird.df$area)
bird.df$temp.sc=scale(bird.df$temp)
bird.df$precip.sc=scale(bird.df$precip)

####
####  Look for Collinearity and Outliers  
####

#pairs(bird.df[,8:11])  # NOTE: There appear to be two outliers
#cor(bird.df[,9:11])

idx.outlier=(1:N)[bird.df$spp==min(bird.df$spp) | bird.df$area==max(bird.df$area)]
bird.sm.df=bird.df[-idx.outlier,]
n=dim(bird.sm.df)[1]

####
####  Create Data Vector and Design Matrix
####

y=bird.sm.df$spp
X=model.matrix(spp~area.sc+temp.sc,data=bird.sm.df)
p=dim(X)[2]

####
####  Fit Poisson Regression Model 
####

source("pois.reg.mcmc.R")
n.mcmc=20000

mu.beta=rep(0,p)
s2.beta=1000^2
beta.tune=.01

out.1=pois.reg.mcmc(y,X,mu.beta,s2.beta,beta.tune,n.mcmc)

layout(matrix(1:3,3,1))
plot(out.1$beta.save[1,],type="l",lty=1)
plot(out.1$beta.save[2,],type="l",lty=1)
plot(out.1$beta.save[3,],type="l",lty=1)

out.1$p.value

####
####  Fit NegBinom Regression Model 
####

source("nb.reg.mcmc.R")
n.mcmc=20000

mu.beta=rep(0,p)
s2.beta=1000^2
beta.tune=.01

mu.logN=log(5)
sig.logN=10

out.2=nb.reg.mcmc(y,X,mu.beta,s2.beta,beta.tune,mu.logN,sig.logN,n.mcmc)

layout(matrix(1:4,4,1))
plot(out.2$beta.save[1,],type="l",lty=1)
plot(out.2$beta.save[2,],type="l",lty=1)
plot(out.2$beta.save[3,],type="l",lty=1)
plot(log(out.2$N.save),type="l",lty=1)

out.2$p.value

####
####  Model Inference Comparison
####

layout(matrix(1:4,2,2))
plot(density(out.1$beta.save[1,]),lwd=2,col=rgb(0,0,0,.5),main="",xlab=bquote(beta[0]))
lines(density(out.2$beta.save[1,]),lwd=2,col=rgb(1,0,0,.5))
curve(dnorm(x,mu.beta[1],sqrt(s2.beta)),lwd=2,col=rgb(0,1,0,.5),add=TRUE)
plot(density(out.1$beta.save[2,]),lwd=2,col=rgb(0,0,0,.5),main="",xlab=bquote(beta[1]))
lines(density(out.2$beta.save[2,]),lwd=2,col=rgb(1,0,0,.5))
curve(dnorm(x,mu.beta[2],sqrt(s2.beta)),lwd=2,col=rgb(0,1,0,.5),add=TRUE)
plot(density(out.1$beta.save[3,]),lwd=2,col=rgb(0,0,0,.5),main="",xlab=bquote(beta[2]))
lines(density(out.2$beta.save[3,]),lwd=2,col=rgb(1,0,0,.5))
curve(dnorm(x,mu.beta[3],sqrt(s2.beta)),lwd=2,col=rgb(0,1,0,.5),add=TRUE)
legend("topleft",lwd=2,col=c(rgb(0,0,0,.5),rgb(1,0,0,.5),rgb(0,1,0,.5)),legend=c("Pois","NB","Prior"))
plot(density(log(out.2$N.save)),lwd=2,col=rgb(1,0,0,.5),main="",xlab="log(N)")
curve(dnorm(x,mu.logN,sig.logN),lwd=2,col=rgb(0,1,0,.5),add=TRUE)



