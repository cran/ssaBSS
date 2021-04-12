\name{rtvvar}
\encoding{latin1}
\alias{rtvvar}
\title{
Simulation of Time Series with Time-varying Variance
}

\description{
Simulating time-varying variance based on TV-VAR model
}

\usage{
rtvvar(n, alpha, beta = 1, simple = FALSE)
}

\arguments{
  \item{n}{Length of the time series}
  \item{alpha}{Parameter \eqn{\alpha} in TV-VAR}
  \item{beta}{Parameter \eqn{\beta} in TV-VAR. Default is 1.}
  \item{simple}{A logical vector indicating whether \eqn{h_t} is considered as its own process, or just \eqn{t/T}. Default is FALSE.}
}

\details{
Time varying variance (TV-VAR) process \eqn{x_t} with parameters \eqn{\alpha} and \eqn{\beta} is of the form
\deqn{x_t = \tilde h_t \epsilon_t,}{x_t = hm_t \epsilon_t,}
where, if \code{simple = FALSE},
\deqn{\tilde h_t^2 = h_t^2 + \alpha x_{t-1}^2,}{hm_t^2 = h_t^2 + \alpha x_(t-1)^2,} where \eqn{\epsilon} are iid \eqn{N(0,1)}, \eqn{x_0 = 0} and \eqn{h_t =  10 - 10 \sin(\beta \pi t/T + \pi/6)  (1 + t/T)}{h_t =  10 - 10 sin(\beta \pi t/T + \pi/6)  (1 + t/T)},

and if \code{simple = TRUE}, \deqn{\tilde h_t = t/T.}{hm_t = t/T.}
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
X <- rtvvar(n, alpha = 0.2, beta = 0.5, simple = FALSE)
plot(X)

}

\keyword{ datagen }
\keyword{ ts }