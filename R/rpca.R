

rpca_alm <- function(M, tol=1e-6) {
  mdim <- dim(M)

  mu <- prod(mdim) / (4.0 * sum(abs(M)))
  lamb <- 1/sqrt(max(mdim))

  L <- matrix(0, dim(M))
  S <- matrix(0, dim(M))
  Y <- matrix(0, dim(M))

  while (!converged(M,L,S,tol)) {
    L <- svd_shrink(M - S - (1/mu) * Y, mu)
    S <- shrink(M - L + (1/mu) * Y, lamb * mu)
    Y <- Y + mu * (M - L - S)
  }

}

shrink <- function(X, tau) {
  aX <- abs(X) - tau
  sign(X) * pmax(aX,0)
}

svd_shrink <- function(X, tau) {
  res <- svd(X)
  res$u %*% (diag(shrink(res$d, tau)) %*% res$v)
}

converged <- function(M,L,S, tol) {
  error = norm(M - L - S, "F") / norm(M, "F")
  message("error =", error)
  error <= tol
}

