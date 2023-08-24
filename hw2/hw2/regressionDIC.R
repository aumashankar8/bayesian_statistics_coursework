norm.DIC <- function(y,X,sigma2,beta.mn,beta.var,n.mcmc,no.print=FALSE){
  library(mvtnorm)
  n <- length(y)
  p <- length(beta.mn)
  
  #Setup Hyperparameters
  n.burn <- round(.1*n.mcmc)
  mu_log <- 0
  s2_log <- 2.25
  Sig.beta <- beta.var*diag(p)
  
  #Initialize outputs
  beta.save <- matrix(0,p,n.mcmc)
  s2.save <- rep(0,n.mcmc)
  Dbar.save <- rep(0,n.mcmc)
  y.pred.mn <- rep(0,n)
  
  #Start
  beta <- solve(t(X)%*%X)%*%t(X)%*%y
  
  ###
  ### MCMC Loop
  ### (COPIED FOR REFERENCE TO SEE WHATS HELPFUL)
  
  for(k in 1:n.mcmc){
    if(!no.print){ if(k%%1000==0) cat(k," ") }
    
    ###
    ### Sample log(s2)
    ###
    # Use samples from Earlier to Sample Beta
    s2 <- sigma2
    
    ###
    ### Sample beta
    ###
    
    tmp.var=solve(t(X)%*%X/s2[k] + solve(Sig.beta))
    tmp.mn=tmp.var%*%(t(X)%*%y/s2[k] + solve(Sig.beta)%*%beta.mn)
    
    beta=as.vector(rmvnorm(1,tmp.mn,tmp.var,method="chol"))
    
    ###
    ### DIC Calculations 
    ###
    
    Dbar.save[k]=-2*sum(dnorm(y,X%*%beta,sqrt(s2[k]),log=TRUE))
    
    ###
    ### Posterior Predictive Calculations 
    ###
    
    if(k > n.burn){
      y.pred=rnorm(n,X%*%beta,sqrt(s2))
      y.pred.mn=y.pred.mn+y.pred/(n.mcmc-n.burn)
    }
    ###
    ### Save Samples
    ###
    
    beta.save[,k]=beta
  }
  
  if(dim(X)[2]==1){
    postbetamn=mean(beta.save[,-(1:n.burn)])
  }
  if(dim(X)[2]>1){
    postbetamn=apply(beta.save[,-(1:n.burn)],1,mean)
  }
  posts2mn=mean(sigma2[-(1:n.burn)])
  Dhat=-2*(sum(dnorm(y,X%*%postbetamn,sqrt(posts2mn),log=TRUE)))
  Dbar=mean(Dbar.save[-(1:n.burn)])
  pD=Dbar-Dhat
  DIC=Dhat+2*pD
  
  if(!no.print){
    cat("Posterior Mean for Beta:","\n")
    print(postbetamn)
    cat("Posterior Mean for s2:","\n")
    print(posts2mn)
    cat("Dhat:",Dhat,"Dbar:",Dbar,"pD:",pD,"DIC:",DIC,"\n")
  }
  
  ###
  ###  Write Output
  ###
  
  list(beta.save=beta.save,s2.save=s2.save,y=y,X=X,n.mcmc=n.mcmc,n=n,r=r,q=q,p=p,Dhat=Dhat,Dbar=Dbar,pD=pD,DIC=DIC,y.pred.mn=y.pred.mn)
  
}
  