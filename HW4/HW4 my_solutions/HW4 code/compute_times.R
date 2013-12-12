#HW 4, Problem 1
#Eugene Shvarts
#
#quick no-frills code to check the execution time. The multiple write-to-file
#routine needs to be fixed -- as it stands, inputting all files into R just
#repeatedly overwrites *PU_time, so I had to do some manual labor there.

library(RCUDA)
source('utility.R')       #loads the function compute_grid

m <- loadModule('rtruncnorm.ptx')
mykernel <- m$rtruncnorm_kernel

NUM_RUNS = 7
MU_VAL = 2
SIG_VAL = 1
LO_VAL = 0
HI_VAL = 1.5
RNG_VAL = c(1L,2L,3L)
NAIVE_VAL = 500L
MAX_VAL = 1000L

#GPU_vals <- rep(0,NUM_RUNS)

for (i in 1:NUM_RUNS) {
  n <- as.integer(10^i)
  x <- rep(0,0,n)
  mu <- rep(MU_VAL,n)
  sigma <- rep(SIG_VAL,n)
  lo <- rep(LO_VAL,n)
  hi <- rep(HI_VAL,n)
  A <- compute_grid(n)
  
  GPU_time <- system.time({
    .cuda(mykernel, out = x, n , mu, sigma, lo, hi, RNG_VAL[1], 
          RNG_VAL[2], RNG_VAL[3], NAIVE_VAL, MAX_VAL, blockDim = A$block_dims, 
          gridDim = A$grid_dims, outputs = 'out')
  })
  fname = paste('GPU_time_',i,'.RData',sep="")
  save(GPU_time, file = fname)
}

source('rtrunc_cpu.R')        #loads the CPU versions of the rtrunc function

#CPU_vals <- rep(0,NUM_RUNS)

for (i in 1:NUM_RUNS) {
  n <- 10^i
  CPU_time <- system.time({
    rtrunc_cpu(n,MU_VAL,SIG_VAL,LO_VAL,HI_VAL,NAIVE_VAL,MAX_VAL)
  })
  fname = paste('CPU_time_',i,'.RData',sep="")
  save(CPU_time, file = fname)
}