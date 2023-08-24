logit.reg.mcmc <- function(y,X,betamn,betavar,beta.tune,n.mcmc){

#
#  Logistic regression with uncorrelated errors
#
#  EXAMPLE USE:
#  ------------
#
#  tmp.out=logit.reg.mcmc(y,X,rep(0,3),100,.1,1000)
#

###
###  Libraries and Subroutines 
###

logit <- function(theta){
  log(theta/(1-theta))
}

logit.inv <- function(z){
  exp(z)/(1+exp(z))
}


###
###  Preliminary Variables
###

n.burn=round(n.mcmc/10)
X=as.matrix(X)
y=as.vector(y)
n=length(y)
p=dim(X)[2]

betasave=matrix(0,p,n.mcmc)
Davgsave=rep(0,n.mcmc)

###
###  Starting Values
###

#beta=matrix(glm(y~X[,-1],family=binomial(link="logit"))$coefficients,p,1)
beta=betamn
theta=logit.inv(X%*%beta)

###
###  Gibbs Loop
###

for(k in 2:n.mcmc){
  if(k%%1000==0) cat(k," ")

  ###
  ### Sample beta
  ###

  beta.star=rnorm(p,beta,beta.tune)
  theta.star=logit.inv(X%*%beta.star)
  mh1=sum(dbinom(y,1,theta.star,log=TRUE))-t(beta.star-betamn)%*%(beta.star-betamn)/betavar/2
  mh2=sum(dbinom(y,1,theta,log=TRUE))-t(beta-betamn)%*%(beta-betamn)/betavar/2
  mh=exp(mh1-mh2)
  if(mh > runif(1)){
    beta=beta.star
    theta=theta.star
  }

  ###
  ### DIC Calculations 
  ###

  Davgsave[k]=-2*(sum(dbinom(y,1,theta,log=TRUE)))

  ###
  ### Save Samples 
  ###

  betasave[,k]=beta

}
cat("\n")

###
###  Calculate DIC 
###

postbetamn=apply(betasave[,-(1:n.burn)],1,mean)
print("Posterior Mean for Beta:")
print(postbetamn)
Dhat=-2*(sum(dbinom(y,1,logit.inv(X%*%postbetamn),log=TRUE)))
Davg=mean(Davgsave[-(1:n.burn)])
pD=Davg-Dhat
DIC=2*Davg-Dhat

cat("Dhat:",Dhat,"Davg:",Davg,"pD:",pD,"DIC:",DIC,"\n")
 
###
###  Write output 
###

list(y=y,X=X,n.mcmc=n.mcmc,betasave=betasave)

}
