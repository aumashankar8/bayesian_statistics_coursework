ar1.mcmc <- function(y,n.mcmc){

#
#
#  Fits AR1 model to stationary data (automatically centers). 
#
#  Example Use:
#  tmp.out=ar1.mcmc(y,1000)
#

###
###  Subroutines 
###

dIG <- function(x,r,q){
  x^(-(q+1))*exp(-1/r/x)/(r^q)/gamma(q)
}

###
###  Center Data
###

T=length(y) 
z=y-mean(y)

###
###  Variables, Priors, Starting Values 
###

s2.save=rep(1,n.mcmc)
alpha.save=rep(0,n.mcmc)

alpha=.1
mu.alpha=0
s2.alpha=1
r=1
q=2

###
###  Gibbs Loop 
###
for(k in 1:n.mcmc){
  if(k%%1000==0) cat(k," ")

  ###
  ###  Sample s2
  ###

  r.tmp=1/(sum((z[-1]-alpha*z[-T])^2)/2 + 1/r)
  q.tmp=(T-1)/2 + q
  s2=1/rgamma(1,q.tmp,,r.tmp)

  ###
  ###  Sample alpha 
  ###

  tmp.var=1/(sum(z[-T]^2)/s2 + 1/s2.alpha)
  tmp.mn=tmp.var*(sum(z[-1]*z[-T])/s2 + mu.alpha/s2.alpha)
  alpha=rnorm(1,tmp.mn,sqrt(tmp.var))

  ###
  ###  Save Samples 
  ###
  
  alpha.save[k]=alpha 
  s2.save[k]=s2

}
cat("\n")
###
###  Make Plots 
###

layout(matrix(1:3,3,1))
plot(z,type="l",lwd=2,ylab="z",xlab="time",main="time series")
abline(h=0,col=8)
plot(density(alpha.save),type="l",lwd=2,xlab="alpha",main="",xlim=c(-1,1))
curve(dnorm(x,mu.alpha,sqrt(s2.alpha)),add=TRUE,col=3,lwd=2)
plot(density(s2.save),type="l",lwd=2,xlab="s2",main="")
curve(dIG(x,r,q),add=TRUE,col=3,lwd=2)

###
###  Write Output
###

list(z=z,s2.save=s2.save,alpha.save=alpha.save,n.mcmc=n.mcmc)

}
