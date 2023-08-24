#Implied Prior?
logit<-function(theta){log(theta)-log(1-theta)}
logit.inv<-function(x){exp(x)/(1+exp(x))}

rnorm(10000, 0, sqrt(2.25))
hist(rnorm(10000, 0, sqrt(2.25)))
hist(logit.inv(rnorm(100000, 0, sqrt(2.25)))) #more in middle
hist(logit.inv(rnorm(100000, 0, sqrt(100)))) #Bathtub shape no good from high sigma^2
