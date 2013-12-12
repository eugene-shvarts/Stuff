#HW 4, Problem 1
#Eugene Shvarts
#
#A code to sample a truncated normal on CPU; I use the same method as for the
#GPU, cf. Robert (2009). See the lecture notes or homework.

"trunc_1side" <- function(x,a_opt,sgn,offset) {   #x ~ U[0,1]
  x <- sgn * offset - log(x) / a_opt               #now an exponential RV
  temp <- exp(-0.5 * (x-a_opt)^2)
  u <- runif(length(x),0,1)
  acceptx <- u <= temp
  return (data.frame(x,acceptx))
}

"trunc_2side" <- function(x,a_opt,mu_m,mu_p) {    #x ~ U[0,1]
  x <- (mu_p-mu_m) * x + mu_m                      #now a rescaled uniform RV
  temp <- exp(0.5 * (a_opt^2 - x^2))
  u <- runif(length(x),0,1)
  acceptx <- u <= temp
  return (data.frame(x,acceptx))
}


"rtrunc_cpu" <- function(n=10000,mu=0,sigma=1,a,b,maxnaive=500,maxtries=1000,verbose=FALSE) {
  ind <- rep(TRUE,n) #ind keeps track of which samples still need to be assigned.
  tind <- which(ind) #supplies references to TRUE elements of ind.
  N <- length(tind)
  samples <- rep(NA,n)
  #coerce all scalars / short vectors:
  mu <- rep_len(mu,N); sigma <- rep_len(sigma,N); a <- rep_len(a,N); b <- rep_len(b,N)
  
  #Let's try the naive method first.
  
  for (i in 1:maxnaive) {               
    x <- rnorm(N,mu[ind],sigma[ind])            #Generate some rnorms
    goodinds <- which(x>= a[ind] & x<= b[ind])  #Keep the useful ones
    samples[tind[goodinds]] <- x[goodinds]
    ind[tind[goodinds]] <- FALSE ; 
    tind <- which(ind)
    N <- length(tind)
    if (N==0) {
      return(samples)
    }
    if(verbose) cat('Naive attempt',i,', N =',N,'\n')
  }
  #Looks like the naive method failed; let's use Robert sampling.
  #Note at this point we are certain that at least one side is truncated.
  
  afin <- tind[which(is.finite(a[ind]))] ; ifelse(length(afin)==0,skipa<-TRUE,skipa<-FALSE)
  bfin <- tind[which(is.finite(b[ind]))] ; ifelse(length(bfin)==0,skipb<-TRUE,skipb<-FALSE)
  bothfin <- intersect(afin,bfin)        ; ifelse(length(bothfin)==0,skipboth<-TRUE,skipboth<-FALSE)
  aexcl <- setdiff(afin,bothfin)         ; ifelse(length(aexcl)==0,skipax<-TRUE,skipax<-FALSE)
  bexcl <- setdiff(bfin,bothfin)         ; ifelse(length(bexcl)==0,skipbx<-TRUE,skipbx<-FALSE) 
  if (length(union(afin,bfin)) != N) 
    cat("Something's wrong -- an untruncated normal wasn't sampled.\n")
# Now we know exactly which side(s) is/are truncated.
  
  mu_m = mu_p <- rep(NA,n)
  mu_m[afin] <- (a[afin] - mu[afin])/sigma[afin]
  mu_p[bfin] <- (b[bfin] - mu[bfin])/sigma[bfin]
  
  sgn <- rep(1,n)
  if (!skipbx)
    sgn[bexcl] <- -1        #sign parameter to switch the one-sided truncation
  
  
  a_opt <- rep(NA,n)       #optimal-alpha parameter
  if(!skipboth)
    a_opt[bothfin] <- apply(cbind(mu_m[bothfin],mu_p[bothfin],0),1,sort)[2,] #choose the middle parameter; 
  a_opt[aexcl] <- 0.5 * (mu_m[aexcl] + sqrt(mu_m[aexcl]^2 + 4))              #see Robert (2009).
  a_opt[bexcl] <- 0.5 * (-mu_p[bexcl] + sqrt(mu_p[bexcl]^2 + 4))
  
  #Necessary parameters generated; now the fun part.
  for (i in 1:maxtries) {

    if(!skipboth) {
        if(verbose) cat('Robert attempt',i,', part 1, N =',N,'\n')
      x <- runif(length(bothfin),0,1)     #First RVs for the two-sided case ...               
      x <- trunc_2side(x,a_opt[bothfin],mu_m[bothfin],mu_p[bothfin])
      goodinds <- bothfin[which(x[,2])]
      goodvals <- x[which(x[,2]),1]     #Keep the useful ones
      samples[goodinds] <- mu[goodinds] + sgn[goodinds] * sigma[goodinds] * goodvals
      # Recall we've been working with standard normals. 
      ind[goodinds] <- FALSE
      tind <- which(ind)
      N <- length(tind)
      if (N==0) return(samples)
      afin <- intersect(afin,tind)
    }
    
    if(!skipa) {
        if(verbose) cat('Robert attempt',i,', part 2, N =',N,'\n')
      x <- runif(length(afin),0,1)        #Then RVs for the right case ...               
      x <- trunc_1side(x,a_opt[afin],sgn[afin],mu_m[afin])
      goodinds <- afin[which(x[,2])]     
      goodvals <- x[which(x[,2]),1]       #Keep the useful ones
      samples[goodinds] <- mu[goodinds] + sgn[goodinds] * sigma[goodinds] * goodvals
      # Recall we've been working with standard normals. 
      ind[goodinds] <- FALSE
      tind <- which(ind)
      N <- length(tind)   
      if (N==0) return(samples)
      bfin <- intersect(bfin,tind)
    }

    if(!skipb){
        if(verbose) cat('Robert attempt',i,', part 3, N =',N,'\n')
      x <- runif(length(bfin),0,1)        #Finally RVs for the left case ...               
      x <- trunc_1side(x,a_opt[bfin],sgn[bfin],mu_p[bfin])
      goodinds <- bfin[which(x[,2])]     
      goodvals <- x[which(x[,2]),1]       #Keep the useful ones
      samples[goodinds] <- mu[goodinds] + sgn[goodinds] * sigma[goodinds] * goodvals
      # Recall we've been working with standard normals. 
      ind[goodinds] <- FALSE
      tind <- which(ind)
      N <- length(tind)   
      if (N==0) return(samples)
      bothfin <- intersect(bothfin,tind)
    }
  }    
  #Oops, Robert sampling failed too. Print a warning message and return NA for missed samples.
  
  cat('Warning:',n,'samples asked for, but only',n-N,'could be obtained. Try increasing maxtries.\n')
  return(samples)
}