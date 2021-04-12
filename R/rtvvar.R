# Simulate time series with time-varying variance
rtvvar <- function(n, alpha, beta = 1, simple = FALSE) {
  x <- h <- hm <- numeric(n) #Zero vectors
  for (i in 2:n) { #The first value stays zero
    if (simple == FALSE) {
    h[i] <- h_fun(i, beta, n)
    hm[i] <- sqrt(h[i]^2 + alpha * x[i-1]^2)
    x[i] <- hm[i] * rnorm(1) 
    } else {
      x[i] <- (i/n) * rnorm(1)       
    }
  }
  as.ts(x)
}

h_fun <- function(t, beta, n) {
  10 - 10 * sin(beta*pi*t/n + pi/6) * (1 + t/n)
}

