#HW 4, Problem 1
#Eugene Shvarts
#
#Code for execution of either the GPU or CPU sampling of rtruncnorm. 
#I did the scalar-to-vector "optimization" in a lazy way; the GPU doesn't see 
#any improvement, but the user can enter scalars or smaller-than-n vectors, 
#non-specified ints, etc. without fear (as long as they anticipate the outcome).
#
#I presume the GPU improvement would involve making each call to mu, sigma, etc.
#the value of idx mod len(var) instead of idx; I don't know whether the GPU can
#interpret this correctly though so I'll leave it to when I have less time
#pressure.

library(RCUDA)
source('utility.R')       #loads the function compute_grid

m <- loadModule('rtruncnorm.ptx') #loads the (pre-compiled) GPU rtrunc function
mykernel <- m$rtruncnorm_kernel   #grab the kernel so we may call it

source('rtrunc_cpu.R')        #loads the CPU versions of the rtrunc function

#Will rep each argument to its needed size.
#No warning provided if the initial size doesn't divide n.
"exec_GPU_rt" <- function(n,mu,sigma,lo,hi,rnga = 1L,rngb = 2L,rngc = 3L, naive = 500L, tries = 1000L) {    
  mu = rep_len(mu,n)
  sigma = rep_len(sigma,n)
  lo = rep_len(lo,n)
  hi = rep_len(hi,n)
  x = rep_len(0.0,n)
  rnga = as.integer(rnga)    #in case user defines these.
  rngb = as.integer(rngb)
  rngc = as.integer(rngc)
  naive = as.integer(naive)
  tries = as.integer(tries)
  n = as.integer(n)
  A <- compute_grid(n)
  
  samples <- .cuda(mykernel, out = x, n , mu, sigma, lo, hi, rnga, rngb, rngc, 
                   naive, tries, blockDim = A$block_dims, gridDim = A$grid_dims, 
                   outputs = 'out')
  return(samples)
}