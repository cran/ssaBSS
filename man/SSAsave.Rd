\name{SSAsave}
\alias{SSAsave}
\alias{SSAsave.default}
\alias{SSAsave.ts}
\alias{SSAsave.xts}
\alias{SSAsave.zoo}
\title{
Identification of Non-stationarity in Variance
}

\description{
SSAsave method for identifying non-stationarity in variance

}

\usage{
SSAsave(X, \ldots)

\method{SSAsave}{default}(X, K, n.cuts = NULL, \ldots)
\method{SSAsave}{ts}(X, \ldots)
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
\deqn{{\bf M} = \sum_{i = 1}^K \frac{T_i}{T}({\bf I} - {\bf S}_{T_i}) ({\bf I} - {\bf S}_{T_i})^T,}{M = sum[(T_i/(T)) (I - S_Ti)*(I - S_Ti)'],}
where \if{html}{\eqn{i = 1, \ldots, K},} \eqn{K} is the number of breakpoints, \eqn{{\bf I}}{I} is an identity matrix, and \eqn{{\bf S}_{T_i}}{S_Ti} is the sample variance of values of \eqn{{\bf Y}}{Y} which belong to a disjoint interval \eqn{T_i}.

The algorithm finds an orthogonal matrix \eqn{{\bf U}}{U} via eigendecomposition
\deqn{{\bf M} = {\bf UDU}^T.}{M = UDU^T.}
The final unmixing matrix is then \eqn{{\bf W} = {\bf U S}^{-1/2}}{W = U S^(-1/2)}. The first \eqn{k} rows of \eqn{{\bf U}}{U} are the eigenvectors corresponding to the non-zero eigenvalues and the rest correspond to the zero eigenvalues. In the same way, the first \eqn{k} rows of \eqn{{\bf W}}{W} project the observed time series to the subspace of components with non-stationary variance, and the last \eqn{p-k} rows to the subspace of components with stationary variance. 


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
  \item{method }{Name of the method ("SSAsave"), to be used in e.g. screeplot.}

}

\references{
Flumian L., Matilainen M., Nordhausen K. and Taskinen S. (2021) \emph{Stationary subspace analysis based on second-order statistics}. Submitted. Available on arXiv: https://arxiv.org/abs/2103.06148
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

z1 <- rtvvar(n, alpha = 0.2, beta = 0.5)
z2 <- rtvvar(n, alpha = 0.1, beta = 1)
z3 <- arima.sim(n, model = list(ma = c(0.72, 0.24))) 
z4 <- arima.sim(n, model = list(ar = c(0.34, 0.27, 0.18)))

Z <- cbind(z1, z2, z3, z4)
X <- as.ts(tcrossprod(Z, A)) # An mts object

res <- SSAsave(X, K = 6)
res$D # Two non-zero eigenvalues
screeplot(res, type = "lines") # This can also be seen in screeplot
ggscreeplot(res, type = "lines") # ggplot version of screeplot

# Plotting the components as an mts object
plot(res) # The first two are nonstationary in variance


}

\keyword{ multivariate }
\keyword{ ts }
\keyword{ methods }