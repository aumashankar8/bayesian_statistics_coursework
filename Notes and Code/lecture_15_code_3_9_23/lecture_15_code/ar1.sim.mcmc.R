ar1.sim.mcmc <- function(T,alpha.true,n.mcmc){

#
#  Simulates then fits AR1 model to simulated data. 
#
#  Example Use:
#  tmp.out=ar1.sim.mcmc(100,.7,1000)
#

###
###  Subroutines 
###

dIG <- function(x,r,q){
  x^(-(q+1))*exp(-1/r/x)/(r^q)/gamma(q)
}

###
###  Simulate Data
###
 
y=rep(.1,T)
for(t in 2:T){
  y[t]=alpha.true*y[t-1]+rnorm(1)
# y[t]=rnorm(1,alpha.true*y[t-1],1)  
}

###
###  Variables, Priors, Starting Values 
###

s2save=rep(1,n.mcmc)
alphasave=rep(0,n.mcmc)

alpha=.1
mualpha=0
s2alpha=1
r=1
q=2
#r=10
#q=.1

###
###  Gibbs Loop 
###

for(k in 1:n.mcmc){
  if(k%%1000==0) cat(k," ")

  ###
  ###  Sample s2
  ###

  rtmp=1/(sum((y[-1]-alpha*y[-T])^2)/2 + 1/r)
  qtmp=(T-1)/2 + q
  s2=1/rgamma(1,qtmp,,rtmp)

  ###
  ###  Sample alpha 
  ###

  tmpvar=1/(sum(y[-T]^2)/s2 + 1/s2alpha)
  tmpmn=tmpvar*(sum(y[-1]*y[-T])/s2 + mualpha/s2alpha)
  alpha=rnorm(1,tmpmn,sqrt(tmpvar))

  ###
  ###  Save Samples 
  ###
  
  alphasave[k]=alpha 
  s2save[k]=s2

}
cat("\n")
###
###  Make Plots 
###

layout(matrix(1:3,3,1))
plot(y,type="l",lwd=2,ylab="y",xlab="time",main="time series")
plot(density(alphasave),type="l",lwd=2,xlab="alpha",main="",xlim=c(-1,1))
curve(dnorm(x,mualpha,sqrt(s2alpha)),add=TRUE,col=3,lwd=2)
abline(v=alpha.true,col=2,lwd=2)
plot(density(s2save),type="l",lwd=2,xlab="s2",main="")
curve(dIG(x,r,q),add=TRUE,col=3,lwd=2)
abline(v=1,col=2,lwd=2)

###
###  Write Output
###

list(y=y,s2save=s2save,alphasave=alphasave,n.mcmc=n.mcmc)

}
