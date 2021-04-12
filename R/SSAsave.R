# Method SSAsave
SSAsave <- function(X, ...) UseMethod("SSAsave")

# SAVE type for variance non-stationarity identification
SSAsave.default <- function(X, K, n.cuts = NULL, ...) {
  n <- nrow(X)
  p <- ncol(X)
  prep <- BSSprep(X)
  Y <- prep$Y 
  if (is.null(n.cuts)) n.cuts <- ceiling(seq(1, n, length = K + 1))
  N.cuts <- n.cuts + c(rep(0, K), 1)
  COV <- array(0, dim = c(p, p, K))
  for (i in 1:K) {
    YK <- Y[N.cuts[i]:(N.cuts[i + 1] - 1), ]
    slicevar <- crossprod(YK)/nrow(YK)
    COV[, , i] <- nrow(YK)*(diag(p) - slicevar)%*%t(diag(p) - slicevar)/n
  }
  M <- apply(COV, 1:2, sum)
  EVD <- eigen(M, symmetric = TRUE)
  W <- crossprod(EVD$vectors, prep$COV.sqrt.i)
  S <- tcrossprod(prep$X.C, W)
  S <- ts(S, names = paste("Series", 1:p))
  RES <- list(W = W, S = S, M = M, K = K, D = EVD$values, MU = prep$MEAN, n.cut = n.cuts, method = "SSAsave")
  class(RES) <- c("ssabss", "bss")
  RES
}

SSAsave.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SSAsave.default(x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

SSAsave.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SSAsave.default(x, ...)
  S <- xts::as.xts(RES$S)
  attr(S, "index") <- attr(X, "index")
  xts::xtsAttributes(S) <- xts::xtsAttributes(X) #attributes additional to zoo
  RES$S <- S
  RES
}

SSAsave.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SSAsave.default(x, ...)
  S <- zoo::as.zoo(RES$S)
  attr(S, "index") <- attr(X, "index")
  RES$S <- S
  RES
}