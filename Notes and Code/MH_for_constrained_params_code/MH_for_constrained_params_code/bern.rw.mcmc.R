bern.rw.mcmc <- function(y,s.tune,n.mcmc){

###
### Set up variables
###
 
n=length(y)
theta.save=rep(0,n.mcmc)
  
###
### Priors and Starting Value 
###
 
alpha=1 
beta=1 
theta=mean(y)
 
###
### MCMC Loop  
###

for(k in 1:n.mcmc){
  if(k%%1000==0){cat(k," ")} 
 
  ###
  ### Sample theta 
  ###
 
  theta.star=rnorm(1,theta,s.tune) 
  if(theta.star>0 & theta.star<1){
    mh.1=sum(dbinom(y,1,theta.star,log=TRUE))+dbeta(theta.star,alpha,beta,log=TRUE)
    mh.2=sum(dbinom(y,1,theta,log=TRUE))+dbeta(theta,alpha,beta,log=TRUE)
    mh=exp(mh.1-mh.2) 
    if(mh>runif(1)){
      theta=theta.star 
    }
  }
  
  ###
  ### Save Sample 
  ###
 
  theta.save[k]=theta 
   
};cat(k," ") 

###
### Write Output 
###
 
list(theta.save=theta.save,n.mcmc=n.mcmc,alpha=alpha,beta=beta) 
  
}