###
###  Exponential Covariance
###

s2=2
phi=1
d.max=10
d=seq(0,d.max,,100)
cov.gram=s2*exp(-d/phi)

plot(d,cov.gram,lwd=2,col=rgb(0,0,0,.5),type="l",main="",xlab="distance (d)",ylab="covariance")

###
###  Simulate GPs 
###

library(mvnfast)
gray.colors.rev <- function (n, start = 0, end = 1, gamma = 1){ 
  gray(seq.int(to= start^gamma, from = end^gamma, length.out = n)^(1/gamma))
}

n.1=40
n.2=40
s.1=seq(0,d.max,,n.1)
s.2=seq(0,d.max,,n.2)
locs=expand.grid(s.1,s.2)
n=dim(locs)[1]
mu=rep(0,n)

D=as.matrix(dist(locs))
Sig.exp=s2*exp(-D/phi)
Sig.gau=.0000000001*diag(n)+s2*exp(-(D^2)/(phi^2))

set.seed(10001)
y.exp=rmvn(1,mu,Sig.exp)
y.gau=rmvn(1,mu,Sig.gau)

layout(matrix(1:2,1,2))
image(matrix(y.exp,n.1,n.2),y=s.2,x=s.1,asp=TRUE,axes=FALSE,xlab="",ylab="",main="Exponential Cov",col=gray.colors.rev(100))
image(matrix(y.gau,n.1,n.2),y=s.2,x=s.1,asp=TRUE,axes=FALSE,xlab="",ylab="",main="Gaussian Cov",col=gray.colors.rev(100))

###
###  Geostatistical Example Data
###

temp.df=read.csv("20171014_Temp.csv")

y=temp.df$Temp
n=length(y)

locs=as.matrix(temp.df[,2:1])
n.lat=75
n.lon=round(n.lat*diff(range(locs[,1]))/diff(range(locs[,2])))
n.pred=n.lat*n.lon
lon.vals=seq(min(locs[,1]),max(locs[,1]),,n.lon)
lat.vals=seq(min(locs[,2]),max(locs[,2]),,n.lat)
locs.pred=as.matrix(expand.grid(lon.vals,lat.vals))
locs.full=rbind(locs,locs.pred)
locs.full.sc=locs.full
locs.full.sc=scale(locs.full)

X.full=cbind(rep(1,n.pred+n),locs.full.sc,locs.full.sc[,1]^2)
X=X.full[1:n,]
X.pred=X.full[-(1:n),]
p=dim(X)[2]

library(maps)
data(us.cities)
CO.cities=us.cities[c(322,247,360,728),]
CO.cities$name=c("Fort Collins","Denver","Grand Junction","Pueblo")

map("state","Colorado",xlim=c(-110,-101),ylim=c(36.5,41.5))
points(locs,pch=20,cex=3*(.25+(y-min(y))/diff(range(y))),col=rgb(0,0,0,.4))
map("state",add=TRUE)
map.axes()
map.cities(CO.cities, country="CO",pch=18,cex=1.5,label=TRUE)
map.cities(CO.cities, country="CO",pch=18,cex=1.5,label=TRUE,capitals=2)

###
###  Fit Geostatistical Model
###

source("norm.geostat.mcmc.R")
n.mcmc=10000
set.seed(1)
mcmc.out=norm.geostat.mcmc(y=y,X=X,locs=locs,beta.mn=rep(0,p),beta.var=1000,n.mcmc=n.mcmc)

layout(matrix(1:3,3,1))
matplot(t(mcmc.out$beta.save),type="l",lty=1)
plot(mcmc.out$s2.save,type="l")
plot(mcmc.out$phi.save,type="l")

dIG <- function(x,r,q){
  x^(-(q+1))*exp(-1/r/x)/(r^q)/gamma(q)
}

layout(matrix(1:6,3,2))
hist(mcmc.out$beta.save[1,],prob=TRUE,breaks=40,xlab=bquote(beta[0]),main="")
curve(dnorm(x,0,sqrt(1000)),lwd=1.5,add=TRUE,lty=2)
hist(mcmc.out$beta.save[2,],prob=TRUE,breaks=40,xlab=bquote(beta[1]),main="")
curve(dnorm(x,0,sqrt(1000)),lwd=1.5,add=TRUE,lty=2)
hist(mcmc.out$beta.save[3,],prob=TRUE,breaks=40,xlab=bquote(beta[2]),main="")
curve(dnorm(x,0,sqrt(1000)),lwd=1.5,add=TRUE,lty=2)
hist(mcmc.out$beta.save[4,],prob=TRUE,breaks=40,xlab=bquote(beta[3]),main="")
curve(dnorm(x,0,sqrt(1000)),lwd=1.5,add=TRUE,lty=2)
hist(mcmc.out$s2.save,prob=TRUE,breaks=40,xlab=bquote(sigma^2),main="")
curve(dIG(x,mcmc.out$r,mcmc.out$q),add=TRUE,lwd=1.5,lty=2)
hist(mcmc.out$phi.save,prob=TRUE,breaks=40,xlab=bquote(phi),main="")
curve(dgamma(x,mcmc.out$gamma.1,mcmc.out$gamma.2),from=0.001,add=TRUE,lwd=1.5,lty=2)

###
###  Bayesian Kriging  
###

library(mvtnorm)
library(fields)
library(maps)

n.mcmc=100
D=rdist.earth(locs,locs)
D.full=rdist.earth(locs.full,locs.full)
D.pred=rdist.earth(locs.pred,locs.pred)
D.uo=D.full[-(1:n),1:n]
n.pred=n.lat*n.lon
y.pred.save=matrix(0,n.pred,n.mcmc)
set.seed(1)
for(k in 1:n.mcmc){
  if((k%%2)==0) cat(k," ")
  s2.tmp=mcmc.out$s2.save[k]
  phi.tmp=mcmc.out$phi.save[k]  
  beta.tmp=mcmc.out$beta.save[,k]
  Sig.inv.tmp=solve(s2.tmp*exp(-D/phi.tmp))
  tmp.mn=X.pred%*%beta.tmp+s2.tmp*exp(-D.uo/phi.tmp)%*%Sig.inv.tmp%*%(y-X%*%beta.tmp)
  tmp.var=s2.tmp*exp(-D.pred/phi.tmp)-s2.tmp*exp(-D.uo/phi.tmp)%*%solve(s2.tmp*exp(-D/phi.tmp))%*%t(s2.tmp*exp(-D.uo/phi.tmp))
  y.pred.save[,k]=as.vector(rmvnorm(1,tmp.mn,tmp.var,method="chol"))
};cat("\n")

y.pred.mn=apply(y.pred.save,1,mean)
y.pred.sd=apply(y.pred.save,1,sd)

gray.colors.rev <- function (n, start = .2, end = 1, gamma = 1){ 
  gray(seq.int(to= start^gamma, from = end^gamma, length.out = n)^(1/gamma))
}

data(us.cities)
CO.cities=us.cities[c(322,247,360,728),]
CO.cities$name=c("Fort Collins","Denver","Grand Junction","Pueblo")

layout(matrix(1:2,1,2))
map("state","Colorado",xlim=c(-110,-101),ylim=c(36.5,41.5))
title("Pred Mean")
image(matrix(y.pred.mn,n.lon,n.lat),y=lat.vals,x=lon.vals,col=gray.colors.rev(100),asp=TRUE,add=TRUE)
points(locs,pch=20,cex=.5)
map("state",add=TRUE)
map.axes()
map.cities(CO.cities, country="CO",pch=18,cex=1.5,label=TRUE)
map.cities(CO.cities, country="CO",pch=18,cex=1.5,label=TRUE,capitals=2)
map("state","Colorado",xlim=c(-110,-101),ylim=c(36.5,41.5))
title("Pred. SD")
image(matrix(y.pred.sd,n.lon,n.lat),y=lat.vals,x=lon.vals,col=gray.colors.rev(100),asp=TRUE,add=TRUE)
points(locs,pch=20,cex=.5)
map("state",add=TRUE)
map.axes()
map.cities(CO.cities, country="CO",pch=18,cex=1.5,label=TRUE)
map.cities(CO.cities, country="CO",pch=18,cex=1.5,label=TRUE,capitals=2)




