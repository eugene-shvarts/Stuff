#HW 0, Problem 8

#Here we simulate forward using a linear regression model, and use bootstrapping
#to estimate the standard error in the computed mean.
library(MASS)

n <- 100
d <- 3          #inherent dimensionality of x , beta
B <- 1000       #number of bootstrap iterations
beta <- c(1.2,.3,-.9)

xrands <- rnorm(n*(d-1),0,1)
dim(xrands) = c(n,d-1)
x <- cbind(1,xrands)

eps <- rnorm(n,0,1)

y <- rep(0,n)
for (i in 1:n) {
  y[i] <- sum(x[i,]*beta)+eps[i]
}
#Perform bootstrapping to calculate the standard error in the vector of ML 
#estimates for beta
b_hat <- matrix(0,d,B)
for (i in 1:B) {
  ind <- sample(n,n,replace=TRUE)
  xsamp <- x[ind,]
  ysamp <- y[ind]
  b_hat[,i] <- ginv(t(xsamp)%*%xsamp)%*%t(xsamp)%*%ysamp  #Compute the ML estimate
#I used the standard linear formula, with general inverse for predictably
#singular matrices (since many rows are duplicates).
}

b_SEs <- rep(0,d)
for (i in 1:d) {
  b_SEs[i] = sd(b_hat[i,])
}

#Compare to the asymptotic results computed by lm
summary(lm(y ~ x))

#My SE's were 0.09771830 0.09598016 0.10184997
#and lm gave  0.10382    0.10742    0.09734

#My estimates were 1.1074850  0.3754895 -1.1272175
#while lm gave     1.10981    0.37103   -1.12557

#Overall a great agreement!
