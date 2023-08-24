###   Homework 5

###   Libraries
library(rjags)
library(vioplot)
library(readr)
library(mvnfast)
library(forcats)
library(ggridges)
library(ggplot2)
library(coda)

###
###   Read the Data
###

setwd("~/UT Austin/Spring 2023/Bayesian Stats/hw5/hw5") #Set WD
oxboys <- read_csv("oxboys.csv")
#View(oxboys)

n <- max(oxboys$Subject)
J <- sum(oxboys$Subject==1)

#plot(oxboys$age,oxboys$height,type="p",ylab="height",xlab="age")

#Setup Data
y <- matrix(1,n,J)
for (i in 1:n){
  y[i,] <- oxboys[which(oxboys$Subject==i),]$height
}

X <- oxboys$age
X.list=vector("list",n)

for(i in 1:n){
  X.list[[i]]=matrix(1,J,2)  
  X.list[[i]][,2]=X[oxboys$Subject==i]
}

X <- matrix(1,n,J)
for (i in 1:n){
  X[i,] <- oxboys[which(oxboys$Subject == i),]$age
}
X <- cbind(1,t(X))

###
###  Setup Models
###

#Without X.List
ox.jags <-"
model{
for(i in 1:n){
      y[i,] ~ dmnorm(X[,c(1,i+1)]%*%beta,tau*I.matrix)
}
    tau ~ dgamma(a,b)
    beta ~ dmnorm(mu.beta, taubeta.matrix)
    mu.beta ~ dmnorm(mu_u, tau_u)

    tau0 ~ dgamma(a0,b0)
    taubeta.matrix[1,1] <- tau0

    tau1 ~ dgamma(a1,b1)
    taubeta.matrix[2,2] <- tau1
    taubeta.matrix[1,2] <- 0
    taubeta.matrix[2,1] <- 0
}
"
mod<-textConnection(ox.jags)

#Attempt with X.List
ox.jags <-"
  model{
    for(i in 1:n){
      y[i,] ~ dmnorm(X.list[[i]]%*%beta,tau*I.matrix)
    }
    tau ~ dgamma(a,b)
    beta ~ dmnorm(mu.beta, taubeta.matrix)
    mu.beta ~ dmnorm(mu_u, tau_u)
    
    tau0 ~ dgamma(a0,b0)
    taubeta.matrix[1,1] <- tau0

    tau1 ~ dgamma(a1,b1)
    taubeta.matrix[2,2] <- tau1
    taubeta.matrix[1,2] <- 0
    taubeta.matrix[2,1] <- 0
}
"
mod<-textConnection(ox.jags)



n.mcmc=100000
n.burn=round(.2*n.mcmc)


#Set Priors
mu_u=rep(0,2)
I.matrix = diag(J)
taubeta.matrix = matrix(0,2,2)
tau_u=1/1000*diag(2) #Precision. 1 / sigma^2
a  = 1
b  = 3
a0 = 1
b0 = 3
a1 = 1
b1 = 3

m.out<-jags.model(mod,data=list('y'=y, 'X.list'=X.list, 'n'=n, 'I.matrix'=I.matrix,'mu_u'=mu_u,'tau_u'=tau_u,'a'=a,'b'=b,'a0'=a0,'b0'=b0,'a1'=a1,'b1'=b1),n.chains=3,n.adapt=0) 

update(m.out,n.burn) # perform burn-in
m.samples=jags.samples(m.out,c('beta','mu.beta','taubeta.matrix[1,1]','taubeta.matrix[2,2]'),n.mcmc) # fit model post-burn-in


#Trace Plot for Betas and Sigma
beta.post.mat=cbind(m.samples$beta[1,,1],m.samples$beta[2,,1])
mu.beta.post.mat=cbind(m.samples$mu.beta[1,,1],m.samples$mu.beta[2,,1])
sig.post.mat=cbind(m.samples$`taubeta.matrix[1,1]`[1,,1],m.samples$`taubeta.matrix[2,2]`[1,,1]) #Need to do 1/taubetamatrix

matplot(beta.post.mat,type="l",lty=1,ylab=bquote(beta))
matplot(sig.post.mat,type="l",lty=1,ylab=bquote(sigma^2))
matplot(mu.beta.post.mat,type="l",lty=1,ylab=bquote(mu))

#Ridgeline Plots



dist.labels=as.factor(rep(as.character(1:n),each=n.mcmc))
slopes.mat=as.matrix(m.samples[[1]][2,,1])
slopes.mat_mu=as.matrix(m.samples$mu.beta[2,,1])
slopes.df=data.frame(x=c(slopes.mat_mu),y=dist.labels)
ggplot(slopes.df, aes(x=x,y=fct_reorder(dist.labels,rep(1:n,each=n.mcmc))))+geom_density_ridges(rel_min_height = 0.005)+labs(y="posterior density",x="parameter value",title="Slope Coefficients")+geom_vline(xintercept=0,size=1,linetype="dashed",col=rgb(0,0,0,.5))



