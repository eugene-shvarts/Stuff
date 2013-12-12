#HW 4, Problem 1
#Eugene Shvarts
#
#OLD VERSION. This code is cute and fast, but it doesn't pair output with input, 
#so it has no hope of working with non-constant vector mu, sigma, etc.
#
#A code to sample a truncated normal on CPU; I use the same method as for the
#GPU, cf. Robert (2009). See the lecture notes or homework.

"trunc_1side" <- function(x,a_opt,sgn,offset) {   #x ~ U[0,1]
  x = sgn * offset - log(x) / a_opt               #now an exponential RV
  temp = exp(-0.5 * (x-a_opt)^2)
  u = runif(length(x),0,1)
  acceptx <- u <= temp
  return (cbind(x,acceptx))
}

"trunc_2side" <- function(x,a_opt,mu_m,mu_p) {    #x ~ U[0,1]
  x = (mu_p-mu_m) * x + mu_m                      #now a rescaled uniform RV
  temp = exp(0.5 * (a_opt - x^2))
  u = runif(length(x),0,1)
  acceptx <- u <= temp
  return (cbind(x,acceptx))
}


"rtrunc_cpu" <- function(n=10000,mu=0,sigma=1,a=-1,b=1,maxnaive=100,maxtries=300,verbose=FALSE) {
  N = n     #N counts the number of samples left to generate.
  samples = NULL
  #Let's try the naive method first.
  
  for (i in 1:maxnaive) {               
    x <- rnorm(2*N,mu,sigma)            #Generate more than we need
    goodvals <- x[x>= a & x<= b]       #Keep the useful ones
    if (length(goodvals) > N) {
      goodvals <- goodvals[1:N]         #But not too many
    }
    samples = c(samples,goodvals)
    N = N-length(goodvals)
    if(verbose) cat('Naive attempt',i,', N =',N,'\n')
    if (N==0) {
      return(samples)
    }
  }  
  #Looks like the naive method failed; let's use Robert sampling.
  #Note at this point we are certain that at least one side is truncated.
  
  sgn = 0                     #Tells us which side is truncated.
  if (is.finite(a)) {
    offset <- (a-mu)/sigma    #this is essentially mu_m
    sgn <- sgn + 1
  }
  
  if (is.finite(b)) {
    offset <- (b-mu)/sigma    #this is essentially mu_p
    sgn <- sgn - 1
  }
  
  if (sgn==0) {               #get parameters needed for both-sided truncation
    mu_m = (a-mu)/sigma
    mu_p = (b-mu)/sigma
  
    if (mu_p<0) {
      a_opt <- mu_p^2
    }
    else if (mu_m<0) {
      a_opt <- mu_m^2
    }
    else {
      a_opt <- 0
    }
    fhandle <- function(x){trunc_2side(x,a_opt,mu_m,mu_p)}
    sgn = 1   #This is a hack; I need it for when I un-standardize the samples.
  }
  else {        #sgn = 1 or -1; get parameters needed for one-sided truncation
    a_opt <- 0.5 * (sgn * offset + sqrt(offset^2 + 4))
    fhandle <- function(x){trunc_1side(x,a_opt,sgn,offset)}
  }
  
  for (i in 1:maxtries) {
    x <- runif(2*N,0,1)                 #We use these to generate the appropriate RVs
    x <- fhandle(x)                     #Recover the sample attempts and whether they were accepted     
    goodvals <- x[x[,2]==TRUE,1]        #Keep the useful ones
    if (length(goodvals) > N) {
      goodvals <- goodvals[1:N]         #But not too many
    }
    goodvals <- mu + sgn*sigma*goodvals #Recall we've been working with standard normals.
    samples = c(samples,goodvals)
    N = N-length(goodvals)
    if(verbose) cat('Robert attempt',i,', N =',N,'\n')
    if (N==0) {
      return(samples)
    }    
  }
  #Oops, Robert sampling failed too. Print a warning message and return NA for missed samples.
  
  cat('Warning:',n,'samples asked for, but only',n-N,'could be obtained. Try increasing maxtries.\n')
  samples = c(samples,rep(NA,n-N))
  return(samples)
}