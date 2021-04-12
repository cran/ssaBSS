\name{ssabss}
\docType{class}
\alias{ssabss}
\alias{plot.ssabss}
\alias{screeplot.ssabss}
\alias{ggscreeplot.ssabss}
\title{
Class: ssabss
}

\description{
Class 'ssabss' (blind source separation in stationary subspace analysis) with methods plot, screeplot (prints a screeplot of an object of class 'ssabss') and ggscreeplot (prints a screeplot of an object of class 'ssabss' using package \code{\link[ggplot2]{ggplot2}}).

The class 'ssabss' also inherits methods from the class 'bss' in package \code{\link[JADE]{JADE}}: for extracting the components (\code{\link[JADE:bss.components]{bss.components}}), for plotting the components (\code{\link[JADE:plot.bss]{plot.bss}}), for printing (\code{\link[JADE:print.bss]{print.bss}}), and for extracting the coefficients (\code{\link[JADE:coef.bss]{coef.bss}}).

}
\usage{
\method{plot}{ssabss}(x, ...)

\method{screeplot}{ssabss}(x, type = c("lines", "barplot"), xlab = "Number of components",
          ylab = NULL, main = paste("Screeplot for", x$method),
          pointsize = 4, breaks = 1:length(x$D), color = "red", \ldots)
                              
\method{ggscreeplot}{ssabss}(x, type = c("lines", "barplot"), xlab = "Number of components",
            ylab = NULL, main = paste("Screeplot for", x$method),
            pointsize = 4, breaks = 1:length(x$D), color = "red", \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class 'ssabss'.}
  \item{type}{Type of screeplot. Choices are "lines" (default) and "barplot".}
  \item{xlab}{Label for x-axis. Default is "Number of components".}
  \item{ylab}{Label for y-axis. Default is "Sum of pseudo eigenvalues" if method is \code{\link{SSAcomb}} and "Eigenvalues" otherwise.}
  \item{main}{Title of the plot. Default is "Screeplot for ...", where ... denotes for the name of the method used.}
  \item{pointsize}{Size of the points in the plot (for \code{type = "lines"} only). Default is 4.}
  \item{breaks}{Breaks and labels for the x-axis. Default is from 1 to the number of series by 1.}
  \item{color}{Color of the line (if \code{type = "lines"}) or bar (if \code{type = "barplot")}. Default is red.}
  \item{\dots}{Further arguments to be passed to or from methods.}
}
\details{
A screeplot can be used to determine the number of interesting components. For \code{\link{SSAcomb}} it plots the sum of pseudo eigenvalues and for other methods it plots the eigenvalues.
}

\author{
Markus Matilainen
}

\seealso{
 \code{\link{ASSA}}, \code{\link{SSAsir}}, \code{\link{SSAsave}}, \code{\link{SSAcor}}, \code{\link{SSAcomb}}, \code{\link[JADE]{JADE}}, \code{\link[ggplot2]{ggplot2}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes }% use one of  RShowDoc("KEYWORDS")
\keyword{ screeplot }