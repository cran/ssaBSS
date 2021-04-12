# Method ASSA
ASSA <- function(X, ...) UseMethod("ASSA")
# Stationarity in mean and variance: ASSA
ASSA.default <- function(X, K,  n.cuts = NULL, ...) {
  n <- nrow(X)
  p <- ncol(X)
  prep <- BSSprep(X)
  Y <- prep$Y 
  if (is.null(n.cuts)) {
    n.cuts <- ceiling(seq(1, n, length = K + 1))
  } else {
    K <- length(n.cuts) - 1
  }
  SStmats <- array(0, dim = c(p, p, K))
  MMtmats <- array(0, dim = c(p, p, K))
  N.cuts <- n.cuts + c(rep(0, K), 1)
  for (i in 1:K) {
    Yint <- Y[N.cuts[i]:(N.cuts[i + 1] - 1), ]
    nint <- nrow(Yint)
    MMtmats[, , i] <- tcrossprod(colMeans(Yint))
    Smati <- crossprod(Yint)/nint
    SStmats[, , i] <- 0.5 * tcrossprod(Smati)
  }
  M <- (apply(MMtmats, 1:2, sum) + apply(SStmats, 1:2, sum))/K - 0.5*diag(p)
  EVD <- eigen(M, symmetric = TRUE)
  W <- crossprod(EVD$vectors, prep$COV.sqrt.i)
  S <- tcrossprod(prep$X.C, W)
  S <- ts(S, names = paste("Series", 1:p))
  RES <- list(W = W, S = S, M = M, K = K, D = EVD$values, MU = prep$MEAN, n.cut = n.cuts, method = "ASSA")
  class(RES) <- c("ssabss", "bss")
  RES
}

ASSA.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- ASSA.default(x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

ASSA.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- ASSA.default(x, ...)
  S <- xts::as.xts(RES$S)
  attr(S, "index") <- attr(X, "index")
  xts::xtsAttributes(S) <- xts::xtsAttributes(X) #attributes additional to zoo
  RES$S <- S
  RES
}
ASSA.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- ASSA.default(x, ...)
  S <- zoo::as.zoo(RES$S)
  attr(S, "index") <- attr(X, "index")
  RES$S <- S
  RES
}