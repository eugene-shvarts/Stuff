#HW 0, Problem 6

#Here we examine the statistical properties of the output of an AR(1) process.

run_ar <- function(n_howlong = 1000, n_howmanytimes = 200, y0 = 0, rho = 0.9) {
  n <- n_howlong
  I <- n_howmanytimes
  rands <- rnorm(n*I,0,1)
  dim(rands) = c(n,I)       #rands will be the noise parameter epsilon_t
  ar_vals = matrix(0,n,I)
  ar_vals[1,] = y0          #each initial value is y0
  
  for (t in 2:n) {
    ar_vals[t,] = ar_vals[t-1,]*rho+rands[t,]
  }
  return(ar_vals)
}

#Simulate the process with rho = 0.9 and n = 1000; plot the results.
n <- 1000
I <- 200
s <- run_ar(n_howmanytimes=1)
plot(1:n,s)

#Repeat the simulation 200 times; store the result in a matrix.
S <- run_ar()

#Compute the mean and variance across each snapshot of time.
tmeans = rowMeans(S)
tvars = rowMeans(S^2)-tmeans^2

#Plot the mean and variance.
split.screen(c(2,1))
screen(1)
plot(1:n,tmeans,main = 'Mean across each snapshot in time')
screen(2)
plot(1:n,tvars, main = 'Variance across each snapshot in time')
close.screen(all = TRUE)

#Compute the mean and variance across each time series. 
smeans = colMeans(S)
svars = colMeans(S^2)-smeans^2

#Plot the mean and variance.
split.screen(c(2,1))
screen(1)
plot(1:I,smeans,main = 'Mean of each time series')
screen(2)
plot(1:I,svars, main = 'Variance of each time series')
close.screen(all = TRUE)

#Justify the above results theoretically. 
#Maybe when I have more time.