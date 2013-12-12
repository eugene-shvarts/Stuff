# HW 4, Problem 2
# Eugene Shvarts
#
# Code for sampling the regression parameters in Bayesian probit model. Performs
# block Gibbs sampling by sampling a multivariate truncated normal (itself with 
# a Gibbs subroutine.) Code adapted from pseudocode featured in this paper:
#   ** http://ba.stat.cmu.edu/journal/2006/vol01/issue01/held.pdf **

library(lattice)
library(coda)
source('rtrunc_cpu.R')    #Need the truncated-normal sampler!
source('exec_rtrunc.R')   #If CUDA option enabled, this contains exec_GPU_rt


"ind_read" <- function(y) {
# Accepts a vector of y values and returns the corresponding truncation intervals
  out = cbind(y,y)
  out[,1] = ifelse(y == 1,0,-Inf)
  out[,2] = ifelse(y == 1,Inf,0)
  return(out)
}


"rtrunc" <- function(n = 1, mu = 0, sigma = 1, y, doCUDA = FALSE) {
# No-hassle function handle for running rtrunc, regardless of how.
  intervals <- ind_read(y)
  
  if (doCUDA)
    samples <- exec_GPU_rt(n,mu,sigma,intervals[,1],intervals[,2])
  else 
    samples <- rtrunc_cpu(n,mu,sigma,intervals[,1],intervals[,2])
  return (samples)
}
  
  
"probit_mcmc" <- function(X, y, v_inv = matrix(0,ncol(X),ncol(X)), 
                          n_burnin = 500, n_iter = 2000, thin_size = 1,
                          doCUDA = FALSE, doMCMC = TRUE) {
#INPUTS:
#   X: design matrix of input values; n x p , so rows are samples.
#   y: column vector of output values; values must be from {0,1}
#   v_inv: prior precsion attached to regression param beta. mean is assumed 0.
#   n_burnin: how many iterations to run the MCMC for before stability assumed.
#   n_iter: how many MCMC samples to actually capture?
#   thin_size: size of the interval between successive entries to n_iter
#   doCUDA: a logical which decides whether we want to run rtrunc via GPU.
###     TOTAL MCMC RUNS WILL BE n_burnin + n_iter * thin_size     ###
  
#OUTPUTS:
#   beta: an mcmc object containing the sampled params post-burnin

  V <- solve(t(X) %*% X + v_inv)   #full proposal covariance
  L <- t(chol(V))                  #lower Cholesky decomposition of V (for easier sampling)
  S <- V %*% t(X)                  #precompute part of the hat matrix
  
  n <- nrow(X)
  p <- ncol(X)
  
  H = rep(0,n)  #the diagonal of the hat matrix X V X'
  W = rep(0,n)  #will be for efficient mvnorm sampling
  Q = rep(0,n)  #will be for efficient mvnorm sampling
  
  for (i in 1:n) {
    H[i] = X[i,] %*% S[,i]
    W[i] = H[i] / (1-H[i])
    Q[i] = W[i] + 1
  }
  
#Time to initialize the latent variable Z, from truncated standard normal.
  Z <- rtrunc(n=n, y=y, doCUDA = doCUDA)

# B denotes the conditional mean. As beta|Z ~ N(B,V) , we move B around the 
# parameter space, sampling with frequency according to thin_size. The Z update
# demands sampling a truncated multivariate norm, but this is integrated into the
# psuedocode algorithm as part of the Z update. 
  
  B <- S %*% Z
  n_tot = n_burnin + n_iter*thin_size
  beta_samples = matrix(0,p,n_iter)
  
  for (i in 1:n_tot) {
    for (j in 1:n) {        #Subroutine Gibbs to sample from multivar norm
      zold <- Z[j]
      m <- X[j,] %*% B
      m <- m - W[j]*(Z[j]-m)
# Now draw the jth element of Z from an appropriate truncated normal.
      Z[j] <- rtrunc(mu = m, sigma = Q[j], y = y, doCUDA = doCUDA)
      B <- B + (Z[j] - zold) * S[,j]
    }                       #end of subroutine
    if (i > n_burnin) {
      if (i-n_burnin %% thin_size) {
        #Draw a new value of beta!
        beta_samples[,round((i-n_burnin)/thin_size)] <- B + L %*% rnorm(p,0,1)
      }
    }
  } #end of MCMC
  if(doMCMC)
    return(mcmc(t(beta_samples), thin = thin_size))
  else
    return(beta_samples)
}