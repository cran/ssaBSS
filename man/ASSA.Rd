\name{ASSA}
\encoding{latin1}
\alias{ASSA}
\alias{ASSA.default}
\alias{ASSA.ts}
\alias{ASSA.xts}
\alias{ASSA.zoo}
\title{
ASSA Method for Non-stationary Identification
}

\description{
ASSA (Analytic Stationary Subspace Analysis) method for identifying non-stationary components of mean and variance.
}

\usage{
ASSA(X, \ldots)

\method{ASSA}{default}(X, K, n.cuts = NULL, \ldots)
\method{ASSA}{ts}(X, \ldots)
}

\arguments{
  \item{X}{A numeric matrix or a multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}}. Missing values are not allowed.}
  \item{K}{Number of intervals the time series is split into.}
  \item{n.cuts}{A K+1 vector of values that correspond to the breaks which are used for splitting the data. Default is intervals of equal length.}
  \item{\ldots}{Further arguments to be passed to or from methods.}
}

\details{

Assume that a \eqn{p}-variate \eqn{{\bf Y}}{Y} with \eqn{T} observations is whitened, i.e. \eqn{{\bf Y}={\bf S}^{-1/2}({\bf X}_t - \frac{1}{T}\sum_{t=1}^T {\bf X}_{t})}{Y = S^(-1/2)*(X_t - (1/T)*sum_t(X_t))}, \if{html}{for \eqn{t = 1, \ldots, T,}}
where \eqn{{\bf S}}{S} is the sample covariance matrix of \eqn{{\bf X}}{X}.

The values of \eqn{{\bf Y}}{Y} are then split into \eqn{K} disjoint intervals \eqn{T_i}. Algorithm first calculates matrix
\deqn{{\bf M} = \frac{1}{T}\sum_{i = 1}^K \left({\bf m}_{T_i} {\bf m}_{T_i}^T + \frac{1}{2} {\bf S}_{T_i} {\bf S}_{T_i}^T\right) - \frac{1}{2} {\bf I},}{M = (1/T)*sum[m_ti m_Ti' + (1/2)*S_Ti S_Ti'] - I/2,}
where \if{html}{\eqn{i = 1, \ldots, K},} \eqn{K} is the number of breakpoints, \eqn{{\bf I}}{I} is an identity matrix, and \eqn{{\bf m}_{T_i}}{m_Ti} is the average of values of \eqn{{\bf Y}}{Y} and \eqn{{\bf S}_{T_i}}{S_Ti} is the sample variance of values of \eqn{{\bf Y}}{Y} which belong to a disjoint interval \eqn{T_i}.

The algorithm finds an orthogonal matrix \eqn{{\bf U}}{U} via eigendecomposition
\deqn{{\bf M} = {\bf UDU}^T.}{M = UDU^T.}
The final unmixing matrix is then \eqn{{\bf W} = {\bf U S}^{-1/2}}{W = U S^(-1/2)}. The first \eqn{k} rows of \eqn{{\bf U}}{U} are the eigenvectors corresponding to the non-zero eigenvalues and the rest correspond to the zero eigenvalues. In the same way, the first \eqn{k} rows of \eqn{{\bf W}}{W} project the observed time series to the subspace of non-stationary components, and the last \eqn{p-k} rows to the subspace of stationary components. 

}

\value{
 A list of class 'ssabss', inheriting from class 'bss', containing the following components:
  \item{W }{The estimated unmixing matrix.}
  \item{S }{The estimated sources as time series object standardized to have mean 0 and unit variances.}
  \item{M }{Used separation matrix.}
  \item{K }{Number of intervals the time series is split into.}
  \item{D }{Eigenvalues of M.}
  \item{MU }{The mean vector of \code{X}.}
  \item{n.cut }{Used K+1 vector of values that correspond to the breaks which are used for splitting the data.}
  \item{method }{Name of the method ("ASSA"), to be used in e.g. screeplot.}
}

\references{
Hara S., Kawahara Y., Washio T. and von \enc{Bünau}{Bunau} P. (2010). \emph{Stationary Subspace Analysis as a Generalized Eigenvalue Problem}, Neural Information Processing. Theory and Algorithms, Part I, pp. 422-429.

}

\author{
Markus Matilainen, Klaus Nordhausen
}

\seealso{
\code{\link[JADE]{JADE}}
}

\examples{

n <- 5000
A <- rorth(4)

z1 <- arima.sim(n, model = list(ar = 0.7)) + rep(c(-1.52, 1.38), c(floor(n*0.5),
        n - floor(n*0.5)))
z2 <- rtvvar(n, alpha = 0.1, beta = 1)
z3 <- arima.sim(n, model = list(ma = c(0.72, 0.24))) 
z4 <- arima.sim(n, model = list(ar = c(0.34, 0.27, 0.18)))

Z <- cbind(z1, z2, z3, z4)
X <- as.ts(tcrossprod(Z, A)) # An mts object

res <- ASSA(X, K = 6)
screeplot(res, type = "lines") # Two non-zero eigenvalues

# Plotting the components as an mts object
plot(res) # The first two are nonstationary

}

\keyword{ multivariate }
\keyword{ ts }
\keyword{ methods }