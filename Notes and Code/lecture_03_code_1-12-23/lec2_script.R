###
###  Monte Carlo Sampling for Continuous Random Variable
###

alpha=2
beta=1
K=10000

p.vec=rbeta(K,alpha,beta)

hist(p.vec,col=8,breaks=30,xlab="p",prob=TRUE,main="")
curve(dbeta(x,alpha,beta),col=3,lwd=2,add=TRUE)

###
###  Monte Carlo Integration 
###

alpha/(alpha+beta) 	# exact mean
mean(p.vec)		# MC approximation 

alpha*beta/(alpha+beta)^2/(alpha+beta+1) 	# exact variance
var(p.vec)					# MC approximation

###
###  Monte Carlo Sampling for Discrete Random Variable
###

K=10000
N=10
p=.25

y.vec=rbinom(K,N,p) # simulate 10000 binomial realizations

plot(table(c(y.vec,seq(0,N,1)))-1,col=rgb(0,0,0,.5),type="h",xlim=c(0,N),ylab="frequency",xlab="y",lwd=3,main="")

###
###  Monte Carlo Summation
###

p*N 		# exact mean:  E(y|p)=p*N=0.25*10=2.5
mean(y.vec) 	# MC summation to get E(y|p) = sum_{y=0}^N y [y|p]  

p*(1-p)*N 	# exact var:  V(y|p)=p*(1-p)*N=0.25*0.75*10=1.875
var(y.vec)	# MC summation to get Var(y|p) = sum_{y=0}^N (y-E(y|p))^2 [y|p]		

