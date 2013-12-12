#HW 3, Problem 1

# Here we are asked to implement root-finding algorithms using the bisection
# technique, and Newton-Raphson.

getroot_bisection <- function(fhandle,interval,kmax=100,FTOL=1e-12,verbose=FALSE) {
# A bisection algorithm for univariate functions. Not easily extensible because
# it relies on the fact that the boundary of a one-dimensional cube is a finite
# set of points.
#
# INPUTS: 
#   fhandle -- a function handle. The function which we are trying to find a root of. 
#   interval -- a vector of length 2. Provides the endpoints of the interval in which we are searching for 
#     roots. Must satisfy the condition that fhandle(interval[1])*fhandle(interval[2]) < 0. 
#   kmax -- a scalar. The maximum number of iterations to run for.
#   FTOL -- a scalar. Acceptable deviance of the function from 0 for convergence.
#   verbose -- a logical. If TRUE, provides debugging and function readout. 
#
# RETURNS:
#   An estimate of a root in the given interval.

# check that the existence of a root is guaranteed.
  if(!fhandle(interval[1])*fhandle(interval[2])<0) {
    cat('Error: There may not be a root in this interval.')
    return()
  }
  
  x1 <- interval[1]
  x2 <- interval[2]
  s1 <- sign(fhandle(x1))
  
  for(k in 1:kmax) {
    x <- mean(c(x1,x2))
    if(verbose) cat('Iteration',k,', x =',x,'\n')
    fx <- fhandle(x)
    if (abs(fx)<FTOL) {
      if(verbose) cat('Converged after',k,'iterations with an error of',abs(fx),'\n')
      return(x)
    }
    if (fx*s1<0) {
      x2 <- x
    } else {
      x1 <- x
    }
  }
  cat("Error: Doesn't seem to have converged.",kmax,'iterations performed.\n')
  if(verbose) cat('The error is',abs(fx),'\n')
  return(x)
}

getroot_newton <- function(fhandle,Jhandle,x0,kmax=100,FTOL=1e-12,XTOL=1e-12,verbose=FALSE) {
  # A pretty standard Newton-Raphson algorithm for univariate functions. I would
  # write a multivariate code but MatLab does it much more easily. 
  #
  # INPUTS: 
  #   fhandle -- a function handle. The function which we are trying to find a root of. 
  #   Jhandle -- a function handle. The derivative associated with fhandle. 
  #   x0 -- a scalar. The initial guess for the algorithm.
  #   kmax -- a scalar. The maximum number of iterations to run for.
  #   FTOL -- a scalar. Acceptable deviance of function from 0 for convergence.
  #   XTOL -- a scalar. Acceptable deviance between guesses for convergence.
  #   verbose -- a logical. If TRUE, provides debugging and function readout. 
  #
  # RETURNS:
  #   An estimate of a root of the function.
  x <- x0
  fx <- fhandle(x)
  
  for(k in 1:kmax) {
    xold <- x
    x = x - fx/Jhandle(x)
    if(verbose) cat('Iteration',k,', x =',x,'\n')
    if(is.infinite(x)) {
      cat('Error: Singular derivative.\n')
      return()
    }
    
    fx <- fhandle(x)
    if(abs(x-xold)<XTOL || abs(fx)<FTOL) {
      if(verbose) cat('Converged after',k,'iterations with an error of',abs(fx),'\n')
      return(x)
    }
  }
  cat("Error: Doesn't seem to have converged.",kmax,'iterations performed.\n')
  if(verbose) cat('The error is',abs(fx),'\n')
  return(x)
}

# Part C provides a function proportional to a desired likelihood, and asks us to estimate the MLE.
# Because the MLE (for this nice, differentiable function) is a stationary point of the derivative, 
# then we use fhandle = the derivative, and Jhandle = the second derivative. The proportionality is
# irrelevant, as the stationary points / roots are unchanged.

# Likelihood prop. to (2+x)^125 * (1-x)^38 * x^34
# W|A tells me that the derivative is then proportional to:
fhandle = (function(x) (x-1)^37 * x^33 * (x+2)^124 * (197*x^2-15*x-68))
# W|A tells me that the second derivative is then proportional to:
Jhandle = (function(x) 15500 * (x-1)^38 * x^34 * (x+2)^123+500 * (x-1)^37 * x^33 * (36*x-17)*(x+2)^124+2*(x-1)^36 * x^32 * (2556*x^2-2414*x+561)*(x+2)^125)

# Also, the likelihood is an odd-degree polynomial and so does not have a global max.
# Reading into the Rao problem, \theta represents a fraction \in [0,1] and, observing the definition
# of \lambda(\theta), that means \lambda lies in the same interval. Seeing that the likelihood vanishes 
# at each endpoint then gives a well-posed problem.

# By observation, the derivative has only two nontrivial zeros; the roots of the quadratic expression. 
# The rule of signs says only one of these is positive, and by noting that the derivative is an even
# polynomial, and experiences an odd number of sign changes before this positive root, the derivative 
# is increasing at this root (and so the function is maximized). The location of this root is then the MLE
# -- we can proceed without worrying about local optimizers.

mleb <- getroot_bisection(fhandle,c(0.01,0.99))
# The result: \lambda = 0.626821497870982 , correct to 7 digits at 23 iterations, and to 12 digits 
# at 53 iterations. The interval was not precisely [0,1] because this violates the conditions
# for the bisection algorithm (as stated, fhandle evaluates to 0 there).

mlen <- getroot_newton(fhandle,Jhandle,0.6)
# The result: \lambda = 0.626821497870982 , correct to 7 digits at 4 iterations, and to 14 digits
# at 5 iterations. Newton-Raphson is very sensitive to choice of initial condition, and having
# prior knowledge of the root from the bisection run was essential. Running Newton with an initial 
# choice of 0.5 instead of 0.6 converged (slowly) to the minimizer \lambda=0 , and starting at 0.8
# converged to \lambda=1. When it works, it leads to significantly faster (quadratic) convergence.