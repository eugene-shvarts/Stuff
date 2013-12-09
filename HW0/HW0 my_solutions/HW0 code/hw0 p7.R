#HW 0, Problem 7

#Here we demonstrate a use of MC and importance sampling for estimating expectation

#With Z ~ norm(0,1), compute E[exp(-Z^2)] via MC integration. 

N <- 10000000
Zn <- rnorm(N,0,1)
estimate <- sum(exp(-Zn^2))/N
cat('E[exp(-Z^2)] is approximately ',estimate,'\n')
cat('E[exp(-Z^2)] computed by exact integration is 1/sqrt(3) = 0.57730...\n')

#With Z ~ TN(0,1;[-2,1]), compute E[Z] via importance sampling.
#Check the main HW PDF for the analytic computation.

Xn <- runif(N,-2,1)
estimate2 <- sum(Xn*dnorm(Xn))/sum(dnorm(Xn))
cat('E[Z] is approximately ',estimate2,'\n')