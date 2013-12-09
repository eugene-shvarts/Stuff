#HW 0, Problem 2

#Generate the x values
x <- runif(10000,0,2*pi)

#Generate the y values
y <- runif(10000,0,1)

#Generate the polar coordinates
u <- y*cos(x)
v <- y*sin(x)

#Scatterplot of the result
plot(u,v)

#Result looks like uniformly distributed points on the unit circle. No surprise; r = sqrt(u^2+v^2) has the distribution of y, which is unif([0,1]).