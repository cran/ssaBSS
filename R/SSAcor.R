# Method SSAcor
SSAcor <- function(X, ...) UseMethod("SSAcor")

# Changes in autocorrelation
SSAcor.default <- function(X, K, n.cuts = NULL, tau = 1, ...) {
  n <- nrow(X)
  p <- ncol(X)
  prep <- BSSprep(X)
  Y <- prep$Y 
  if (is.null(n.cuts)) n.cuts <- ceiling(seq(1, n, length = K + 1))
  N.cuts <- n.cuts + c(rep(0, K), 1)
  S <- array(0, dim = c(p, p, K))
  Sall <- actau(Y, tau = tau)
  for (i in 1:K) {
    slice <- Y[N.cuts[i]:(N.cuts[i + 1] - 1), ]
    sli.length <- nrow(slice)
    S[, , i] <- sli.length*(Sall - actau(slice, tau = tau))%*%t(Sall - actau(slice, tau = tau))/n
  }
  M <- apply(S, 1:2, sum)
  EVD <- eigen(M, symmetric = TRUE)
  W <- crossprod(EVD$vectors, prep$COV.sqrt.i)
  S <- tcrossprod(prep$X.C, W)
  S <- ts(S, names = paste("Series", 1:p))
  RES <- list(W = W, S = S, M = M, K = K, D = EVD$values, MU = prep$MEAN, n.cut = n.cuts, k = tau, method = "SSAcor")
  class(RES) <- c("ssabss", "bss")
  RES
}

SSAcor.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SSAcor.default(x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

SSAcor.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SSAcor.default(x, ...)
  S <- xts::as.xts(RES$S)
  attr(S, "index") <- attr(X, "index")
  xts::xtsAttributes(S) <- xts::xtsAttributes(X) #attributes additional to zoo
  RES$S <- S
  RES
}

SSAcor.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SSAcor.default(x, ...)
  S <- zoo::as.zoo(RES$S)
  attr(S, "index") <- attr(X, "index")
  RES$S <- S
  RES
}
#auxilliary function
actau <- function(Y, tau) {
  n <- nrow(Y)
  Yt <- Y[1:(n - tau), ]
  Yti <- Y[(1 + tau):n, ]
  crossprod(Yt, Yti)/nrow(Yt)
}
