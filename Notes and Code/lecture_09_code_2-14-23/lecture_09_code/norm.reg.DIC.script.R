####
####  Read in Data
####

med.df=read.csv("MedicalData.csv",header=TRUE)
head(med.df)

####
####  Create Design Matrices
####

y=med.df$BP
n=length(y)

L=4
X.list=vector("list",L)
X.list[[1]]=model.matrix(BP~BMI+Glucose,data=med.df)
X.list[[2]]=model.matrix(BP~BMI,data=med.df)
X.list[[3]]=model.matrix(BP~Glucose,data=med.df)
X.list[[4]]=model.matrix(BP~1,data=med.df)

####
####  Fit Models and Calculate DIC 
####

source("norm.reg.mcmc.R")
beta.mean.list=vector("list",L)
beta.var.list=vector("list",L)
for(l in 1:L){
  beta.mean.list[[l]]=rep(0,dim(X.list[[l]])[2])
  beta.var.list[[l]]=1000
}
mcmc.out.list=vector("list",L)
DIC.vec=rep(0,L)
for(l in 1:L){
  mcmc.out.list[[l]]=norm.reg.mcmc(y,X.list[[l]],beta.mean.list[[l]],beta.var.list[[l]],50000)
  DIC.vec[l]=mcmc.out.list[[l]]$DIC
}

DIC.vec
plot(1:L,DIC.vec,type="o",lwd=2,ylab="DIC",xlab="Model")

####
####  Fit Models and Calculate LOO CV Deviance Score 
####  (This will take several minutes/hours to run!)
####

source("norm.reg.mcmc.R")
beta.mean.list=vector("list",L)
beta.var.list=vector("list",L)
for(l in 1:L){
  beta.mean.list[[l]]=rep(0,dim(X.list[[l]])[2])
  beta.var.list[[l]]=10
}
mcmc.out.list=vector("list",L)
score.vec=rep(0,L)
for(i in 1:n){
  cat(i," ")
  for(l in 1:L){
    cat(l," ")
    mcmc.out.list[[l]]=norm.reg.mcmc(y[-i],X.list[[l]][-i,],beta.mean.list[[l]],beta.var.list[[l]],10000,TRUE)
    score.vec[l]=score.vec[l]-2*mean(dnorm(y[i],t(X.list[[l]][i,])%*%mcmc.out.list[[l]]$beta.save,sqrt(mcmc.out.list[[l]]$s2.save),log=TRUE))
  }; cat("\n")
}

score.vec
plot(1:L,score.vec,type="o",lwd=2,ylab="Score",xlab="Model")


