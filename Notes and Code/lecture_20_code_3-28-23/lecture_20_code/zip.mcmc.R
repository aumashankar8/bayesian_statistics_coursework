zip.mcmc <- function(y,n.mcmc){

#
#  ZIP model with homogeneous lambda
#

###
###  Set up variables
###

n=length(y)

lam.save=rep(0,n.mcmc)
p.save=rep(0,n.mcmc)
z.mean=rep(0,n)

n.burn=round(.2*n.mcmc)

###
###  Priors 
###

alpha.p=1
beta.p=1

alpha.lam=.01
beta.lam=.01

###
###  Starting values 
###

lam=mean(y)
z=rep(0,n)
z[y>0]=1

###
###  Begin MCMC loop 
###

for(k in 1:n.mcmc){
  if(k%%1000==0){cat(k," ")}

  ###
  ###  Sample p
  ###

  p=rbeta(1,sum(z)+alpha.p,sum(1-z)+beta.p)

  ###
  ###  Sample lambda 
  ###
  
  lam=rgamma(1,sum(y[z==1])+alpha.lam,sum(z)+beta.lam)

  ###
  ###  Sample z 
  ###

  p.tmp=p*exp(-lam)/(p*exp(-lam)+1-p)
  z[y==0]=rbinom(sum(y==0),1,p.tmp)

  ###
  ###  Save Samples
  ###

  p.save[k]=p
  lam.save[k]=lam
  if(k>n.burn){
    z.mean=z.mean+z/(n.mcmc-n.burn)
  }

};cat("\n")

###
### Write Output 
###

list(p.save=p.save,lam.save=lam.save,z.mean=z.mean,n.mcmc=n.mcmc)

}
