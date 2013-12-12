// Eugene Shvarts
// STA 250 | Fall 2013

//Accept-reject sampler for the truncated normal distribution (either end finite or infinite).
//Attempts naive rejection sampling until a specified number of failures, then switches to
//the method detailed in Robert (2009) (see the Lecture notes or the assignment).
//In the unlikely event that the Robert method fails, will return NA after a (different)
//specified number of failures. 

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <math_constants.h>


extern "C"
{
__global__ void 
rtruncnorm_kernel(float *vals, int n,             //vals is the input and output
                  float *mu, float *sigma,        //these are the distribution parameters
                  float *lo, float *hi,           //these are the truncation parameters
                  int rng_a, int rng_b, int rng_c,//these are the RNG seeds
                  int maxnaive,                   //after maxnaive attempts, switch methods
                  int maxtries)                   //after maxtries attempts, cancel with error
{
    // Usual block/thread indexing...
    int myblock = blockIdx.x + blockIdx.y * gridDim.x;
    int blocksize = blockDim.x * blockDim.y * blockDim.z;
    int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
    int idx = myblock * blocksize + subthread;

    if (idx < n) {
      // Setup the RNG:
      curandState rng;
      curand_init(rng_a+idx*rng_b, rng_c, 0, &rng);

      // Sample:
      // First try naive rejection method, until maxnaive failures. 
      for (int i = 1; i < maxnaive; i = i+1) {
        vals[idx] = mu[idx] + sigma[idx] * curand_normal(&rng);
        if (vals[idx] > lo[idx] && vals[idx] <= hi[idx]) {
          return;
        }
      }
      // If we made it this far without a return statement, we've had maxnaive failures.
      // So, it's time to try the more sophisticated truncated sampling method cf. Robert.
      //  First we decide between two-sided and one-sided truncation. Note that there must be
      //  some truncation; otherwise it is impossible to reach this point.
      int which_side = 0;
      float offset = 0;
      
      if (isfinite(lo[idx])) {
        offset = (lo[idx]-mu[idx])/sigma[idx];
        which_side = which_side + 1;
      }
      
      if (isfinite(hi[idx])) {
        offset = (hi[idx]-mu[idx])/sigma[idx];
        which_side = which_side - 1;
      }
      //When which_side = 1, truncation is unbounded to the right. 
      //When which_side = -1, unbounded to the left.
      //When which_side = 0, bounded on both sides.
      float a_opt = 0;  //represents the optimal alpha for the exponential distribution
      float temp = 0;   //from sampling the exponential distribution by inverse CDF
      float u = 0;      //generated to compare against the inverse CDF
      
      if (which_side != 0) {  //the single-truncated case
        a_opt = 0.5 * (which_side * offset + sqrtf(powf(offset,2) + 4));
        for (int i = 1; i<maxtries; i = i+1) {
          temp = curand_uniform(&rng);
          vals[idx] = which_side * offset - logf(temp) / a_opt; //the exponential sample
          temp = expf(-0.5 * powf(vals[idx]-a_opt, 2)); //using overwriting to conserve space
          u = curand_uniform(&rng); 
          if (u <= temp) {
            //we used a standardized normal this entire time
            vals[idx] = mu[idx] + which_side*sigma[idx]*vals[idx];
            return;
          }
        } // end maxtries loop
      } // end single-truncated case
      else {  //the both-truncated-sides case
        float mu_m = (lo[idx]-mu[idx])/sigma[idx];
        float mu_p = (hi[idx]-mu[idx])/sigma[idx];
        if (mu_p < 0) {
          a_opt = powf(mu_p,2);
        }
        else if (mu_m > 0) {
          a_opt = powf(mu_m,2);
        }
        for (int i = 1; i<maxtries; i = i+1) {
          vals[idx] = (mu_p-mu_m) * curand_uniform(&rng) + mu_m; //uniform on the appropriate interval
          temp = expf(0.5 * (a_opt - powf(vals[idx],2)));  //choose the right branch of the function
          u = curand_uniform(&rng);
          if (u <= temp) {
            //we used a standardized normal this entire time
            vals[idx] = mu[idx] + sigma[idx]*vals[idx];
            return;
          }
        } // end maxtries loop
      } // end both-truncated case
      // If the code reaches this point, maxtries has been exhausted with no luck. 
    } // end if idx < n
    vals[idx] = CUDART_NAN_F; //so, return NA
    return;
} // end rtruncnorm_kernel

} // END extern "C"

