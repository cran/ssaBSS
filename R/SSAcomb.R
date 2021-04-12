# Method SSAcomb
SSAcomb <- function(X, ...) UseMethod("SSAcomb")

# Combination of MCOR, MSIR and MSAVE: Joint Diagonalization
SSAcomb.default <- function(X, K, n.cuts = NULL, tau = 1, eps = 1e-6, maxiter = 2000, ...) {
  n <- nrow(X)
  p <- ncol(X)
  prep <- BSSprep(X)
  Y <- prep$Y 
  if (is.null(n.cuts)) n.cuts <- ceiling(seq(1, n, length = K + 1))
  N.cuts <- n.cuts + c(rep(0, K), 1)
  
  R <- array(0, dim = c(p, p, 3))
  R[, , 1] <- SSAsir(X, K, n.cuts)$M
  R[, , 2] <- SSAsave(X, K, n.cuts)$M
  R[, , 3] <- SSAcor(X, K, n.cuts, tau)$M
  JD <- frjd(R, eps = eps, maxiter = maxiter)
  D <- JD$D
  sumdg <- diag(apply(D, 1:2, sum))
  ord <- order(sumdg, decreasing = TRUE)
  P <- diag(p)
  P <- P[ord, ]
  DTable <- matrix(0, ncol = p, nrow = 3)
  for (j in 1:3) {
    D[ , , j] <- P %*% tcrossprod(D[ , , j], P) #Diagonal elements are now in order
    DTable[j, ] <- diag(D[ , , j])
  }
  
  rownames(DTable) <- c("Msir", "Msave", "Mcor")
  V <- JD$V[, ord]
  W <- crossprod(V, prep$COV.sqrt.i)
  S <- tcrossprod(prep$X.C, W)
  S <- ts(S, names = paste("Series", 1:p))
  RES <- list(W = W, S = S, R = R, K = K, D = colSums(DTable), DTable = DTable, MU = prep$MEAN, n.cut = n.cuts, k = tau, method = "SSAcomb")
  class(RES) <- c("ssabss", "bss")
  RES
}

SSAcomb.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SSAcomb.default(x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

SSAcomb.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SSAcomb.default(x, ...)
  S <- xts::as.xts(RES$S)
  attr(S, "index") <- attr(X, "index")
  xts::xtsAttributes(S) <- xts::xtsAttributes(X) #attributes additional to zoo
  RES$S <- S
  RES
}

SSAcomb.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SSAcomb.default(x, ...)
  S <- zoo::as.zoo(RES$S)
  attr(S, "index") <- attr(X, "index")
  RES$S <- S
  RES
}

`plot.ssabss` <- function(x, ...) {
  S <- x$S
  if(ncol(S) <= 2) {
    plot(S, ...)
  } else {
    if (any(class(S) %in% c("mts", "xts", "zoo"))) {
      plot(S, ...)
    } else {
      pairs(S, ...)
    }
  }
}
