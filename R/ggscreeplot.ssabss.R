#ggplot2 style screeplot for class ssabss
ggscreeplot.ssabss <- function(x, type = c("lines", "barplot"), xlab = "Number of components", ylab = NULL, main = paste("Screeplot for", x$method),
                                  pointsize = 4, breaks = 1:length(x$D), color = "red", ...) {
  if (is.null(ylab)) {
    if (x$method == "SSAcomb") {
      ylab <- "Sum of pseudo eigenvalues"
    } else {
      ylab <- "Eigenvalues"
    }
  }
  D <- x$D
  idx <- 1:length(D)
  DF <- data.frame(D = D, idx = idx)
  if(type == "lines") {
    ggplot(DF, aes(x = idx, y = D)) + geom_line(color = color) + geom_point(size = pointsize) +
      xlab(xlab) + ylab(ylab) +
      ggtitle(main) +
      theme_bw() + scale_x_continuous(breaks = breaks, labels = breaks)
  } else { #"barplot"
    ggplot(DF, aes(x = idx, y = D)) + geom_bar(stat = "identity", fill = color) +
      xlab(xlab) + ylab(ylab) +
      ggtitle(main) + 
      theme_bw() + scale_x_continuous(breaks = breaks, labels = breaks)
  }
}
