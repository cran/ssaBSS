# Simulate time series with time-varying autovariance - tv-AR1
rtvAR1 <- function(n, sigma = 0.93) {
  x <- numeric(n) #Zero vectors
  for (i in 2:n) {#The first value stays zero
    x[i] <- a_fun(i, n) * x[i-1] + rnorm(1, 0, sigma) 
  }
  as.ts(x)
}

a_fun <- function(t, n) {
  0.5*cos(2*pi*t/n)
}

