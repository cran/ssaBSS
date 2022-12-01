\name{ssaBSS-package}
\encoding{latin1}
\alias{ssaBSS-package}
\docType{package}
\title{
Stationary Subspace Analysis}

\description{
Stationary subspace analysis (SSA) is a blind source separation (BSS) variant where stationary components are separated from non-stationary components. Several SSA methods for multivariate time series are provided here (Flumian et al. (2021); Hara et al. (2010) <doi:10.1007/978-3-642-17537-4_52>) along with functions to simulate time series with time-varying variance and autocovariance (Patilea and Raïssi(2014) <doi:10.1080/01621459.2014.884504>).
}

\details{
\tabular{ll}{
Package: \tab ssaBSS\cr
Type: \tab Package\cr
Version: \tab 0.1.1\cr
Date: \tab 2022-12-01\cr
License: \tab GPL (>= 2)
}
This package contains functions for identifying different types of nonstationarity
\itemize{
\item{\code{\link{SSAsir}}} { -- SIR type function for mean non-stationarity identification}
\item{\code{\link{SSAsave}}} { -- SAVE type function for variance non-stationarity identification}
\item{\code{\link{SSAcor}}} { -- Function for identifying changes in autocorrelation}
\item{\code{\link{ASSA}}} { -- ASSA: Analytic SSA for identification of nonstationarity in mean and variance.}
\item{\code{\link{SSAcomb}}} { -- Combination of \code{\link{SSAsir}}, \code{\link{SSAsave}}, and \code{\link{SSAcor}} using joint diagonalization}
}

The package also contains function \code{\link{rtvvar}} to simulate a time series with time-varying variance (TV-VAR), and function \code{\link{rtvAR1}} to simulate a time series with time-varying autocovariance (TV-AR1).
}


\author{
Markus Matilainen, Léa Flumian, Klaus Nordhausen, Sara Taskinen

Maintainer: Markus Matilainen <markus.matilainen@outlook.com> 
}

\references{
Flumian L., Matilainen M., Nordhausen K. and Taskinen S. (2021) \emph{Stationary subspace analysis based on second-order statistics}. Submitted. Available on arXiv: https://arxiv.org/abs/2103.06148

Hara S., Kawahara Y., Washio T. and von \enc{Bünau}{Bunau} P. (2010). \emph{Stationary Subspace Analysis as a Generalized Eigenvalue Problem}, Neural Information Processing. Theory and Algorithms, Part I, pp. 422-429.

Patilea V. and \enc{Raïssi}{Raissi} H. (2014) \emph{Testing Second-Order Dynamics for Autoregressive Processes in Presence of Time-Varying Variance}, Journal of the American Statistical Association, 109 (507), 1099-1111.

}

\keyword{package}
\keyword{multivariate}
\keyword{ts}