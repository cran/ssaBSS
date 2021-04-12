# screeplot for class ssabss
screeplot.ssabss <- function (x, type = c("lines", "barplot"), xlab = "Number of components",
                              ylab = NULL, main = paste("Screeplot for", x$method), 
                              pointsize = 4, breaks = 1:length(x$D), color = "red", ...) {

    if (is.null(ylab)) {
        if (x$method == "SSAcomb") {
            ylab <- "Sum of pseudo eigenvalues"
        } else {
            ylab <- "Eigenvalues"
        }
    }
    
    if (type == "lines") {
        idx <- 1:length(x$D)
        plot(idx, x$D, type = "l", ylab = ylab, xlab = xlab,
             main = main, xaxt = "n", col = color, ...)
        points(x$D, lwd = pointsize)
        axis(1, at = breaks, labels = TRUE)
    } else {
        labels <- rep("", length(x$D))
        labels[breaks] <- breaks
        barplot(x$D, ylab = ylab, xlab = xlab, main = main, col = color, names.arg = labels, ...)
    }
}

