\name{rtvAR1}
\encoding{latin1}
\alias{rtvAR1}
\title{
Simulation of Time Series with Time-varying Autocovariance
}

\description{
Simulating time-varying variance based on TV-AR1 model
}

\usage{
rtvAR1(n, sigma = 0.93)
}

\arguments{
  \item{n}{Length of the time series}
  \item{sigma}{Parameter \eqn{\sigma^2} in TV-AR1, i.e. the variance. Default is 0.93.}
}

\details{
Time varying autoregressive processes of order 1 (TV-AR1) is
\deqn{x_t = a_t x_{t-1} + \epsilon_t,}{x_t = a_t x_(t-1) + \epsilon_t,}
with \eqn{x_0=0}, \eqn{\epsilon_t} is iid \eqn{N(0, \sigma^2)} and \eqn{a_t = 0.5\cos(2\pi t/T)}.
}

\value{
The simulated series as a \code{\link{ts}} object.
}

\references{
Patilea V. and \enc{Raïssi}{Raissi} H. (2014) \emph{Testing Second-Order Dynamics for Autoregressive Processes in Presence of Time-Varying Variance}, Journal of the American Statistical Association, 109 (507), 1099-1111.

}

\author{
Sara Taskinen, Markus Matilainen
}

\examples{

n <- 5000
X <- rtvAR1(n, sigma = 0.93)
plot(X)

}

\keyword{ datagen }
\keyword{ ts }