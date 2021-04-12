\name{SSAcomb}
\alias{SSAcomb}
\alias{SSAcomb.default}
\alias{SSAcomb.ts}
\alias{SSAcomb.xts}
\alias{SSAcomb.zoo}
\title{
Combination Main SSA Methods
}

\description{
SSAcomb method for identification for non-stationarity in mean, variance and covariance structure.
}

\usage{
SSAcomb(X, \ldots)

\method{SSAcomb}{default}(X, K, n.cuts = NULL, tau = 1, eps = 1e-6, maxiter = 2000, \ldots)
\method{SSAcomb}{ts}(X, \ldots)
}

\arguments{
  \item{X}{A numeric matrix or a multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}}. Missing values are not allowed.}
  \item{K}{Number of intervals the time series is split into.}
  \item{n.cuts}{A K+1 vector of values that correspond to the breaks which are used for splitting the data. Default is intervals of equal length.}
  \item{tau}{The lag as a scalar. Default is 1.}
  \item{eps}{Convergence tolerance.}
  \item{maxiter}{The maximum number of iterations.}
  \item{\ldots}{Further arguments to be passed to or from methods.}
  }

\details{

Assume that a \eqn{p}-variate \eqn{{\bf Y}}{Y} with \eqn{T} observations is whitened, i.e. \eqn{{\bf Y}={\bf S}^{-1/2}({\bf X}_t - \frac{1}{T}\sum_{t=1}^T {\bf X}_{t})}{Y = S^(-1/2)*(X_t - (1/T)*sum_t(X_t))}, \if{html}{for \eqn{t = 1, \ldots, T},}
where \eqn{{\bf S}}{S} is the sample covariance matrix of \eqn{{\bf X}}{X}. 

The values of \eqn{{\bf Y}}{Y} are then split into \eqn{K} disjoint intervals \eqn{T_i}. For a chosen \eqn{\tau}{tau}, algorithm first calculates the \eqn{{\bf M}}{M} matrices from SSAsir (matrix \eqn{{\bf M}_1}{M_1}), SSAsave (matrix \eqn{{\bf M}_2}{M_2}) and SSAcor (matrix \eqn{{\bf M}_3}{M_3}).

The algorithm finds an orthogonal matrix \eqn{{\bf U}}{U} by maximizing 
\deqn{\sum_{i = 1}^3 ||\textrm{diag}({\bf}{\bf U}{\bf M}_i {\bf U}')||^2.}{sum(||diag(U M_i U')||^2),}
\if{html}{where \eqn{i = 1, \ldots, 3}.} The final unmixing matrix is then \eqn{{\bf W} = {\bf US}^{-1/2}}{W = U S^(-1/2)}.

Then the pseudo eigenvalues \eqn{{\bf D}_i = \textrm{diag}({\bf}{\bf U}{\bf M}_i {\bf U}') = \textrm{diag}(d_{i,1}, \ldots, d_{i,p})}{D_i = diag(U M_i U') = (d_i1, ..., d_ip)} are obtained and the value of \eqn{d_{i,j}}{d_ij} tells if the \eqn{j}th component is nonstationary with respect to \eqn{{\bf M}_i}{M_i}.

}

\value{
 A list of class 'ssabss', inheriting from class 'bss', containing the following components:
  \item{W }{The estimated unmixing matrix.}
  \item{S }{The estimated sources as time series object standardized to have mean 0 and unit variances.}
  \item{R }{Used M-matrices as an array.}
  \item{K }{Number of intervals the time series is split into.}
  \item{D }{The sums of pseudo eigenvalues.}
  \item{DTable }{The peudo eigenvalues of size 3*p to see which type of nonstationarity there exists in each component.}
  \item{MU }{The mean vector of \code{X}.}
  \item{n.cut }{Used K+1 vector of values that correspond to the breaks which are used for splitting the data.}
  \item{k }{The used lag.}
  \item{method }{Name of the method ("SSAcomb"), to be used in e.g. screeplot.}
}

\references{
Flumian L., Matilainen M., Nordhausen K. and Taskinen S. (2021) \emph{Stationary subspace analysis based on second-order statistics}. Submitted. Available on arXiv: https://arxiv.org/abs/2103.06148
}

\author{
Markus Matilainen, Klaus Nordhausen
}

\seealso{
\code{\link[JADE]{JADE}} \code{\link[JADE]{frjd}}

}

\examples{

n <- 10000
A <- rorth(6)

z1 <- arima.sim(n, model = list(ar = 0.7)) + rep(c(-1.52, 1.38), 
        c(floor(n*0.5), n - floor(n*0.5)))
z2 <- rtvAR1(n)
z3 <- rtvvar(n, alpha = 0.2, beta = 0.5)
z4 <- arima.sim(n, model = list(ma = c(0.72, 0.24), ar = c(0.14, 0.45)))
z5 <- arima.sim(n, model = list(ma = c(0.34))) 
z6 <- arima.sim(n, model = list(ma = c(0.72, 0.15)))

Z <- cbind(z1, z2, z3, z4, z5, z6)
library(xts)
X <- tcrossprod(Z, A)
X <- xts(X, order.by = as.Date(1:n)) # An xts object

res <- SSAcomb(X, K = 12, tau = 1)
ggscreeplot(res, type = "lines") # Three non-zero eigenvalues
res$DTable # Components have different kind of nonstationarities

# Plotting the components as an xts object
plot(res, multi.panel = TRUE) # The first three are nonstationary


}

\keyword{ multivariate }
\keyword{ ts }
\keyword{ methods }