optshrink <- function(Xtil, r) {

  svdres <- svd(Xtil)
  stil <- svdres$d
  Stil <- diag(svdres$d)

  sigmahats <- stil[1:r]
  Xnoise_est <- Stil[(r+1):nrow(Stil), (r+1):ncol(Stil)]
  theta_hats <- numeric(r)
  Dpz <- numeric(r)
  wopt_hats <- numeric(r)

  for (idx in 1:r) {
    theta_hats[idx] <- sqrt(1/estimateDz(Xnoise_est,sigmahats[idx]))
    Dpz[idx] = estimateDpz(Xnoise_est,sigmahats[idx])
    wopt_hats[idx] = -2/(theta_hats[idx]^2*Dpz[idx])
  }

  Shat = svdres$u[,1:r] %*% diag(wopt_hats) %*% t(svdres$v[,1:r])

  relmse_hat <- 1 - sum(wopt_hats^2)/sum(theta_hats^2);
  mse_hat <- relmse_hat*sum(theta_hats^2)
}

estimateDz <- function(X, z) {
  n <- nrow(X)
  m <- ncol(X)

  In <- diag(n)
  Im <- diag(m)


  z2IXXt = z^2 * In - X %*% t(X)
  z2IXtX = z^2 * Im - t(X) %*% X
  invz2XtX = solve(z2IXtX)
  invz2XXt = solve(z2IXXt)

  D1z = 1/n * sum(diag(z * invz2XXt))
  D2z = 1/m * sum(diag(z * invz2XtX))

  D1z * D2z

}

estimateDpz <- function(X, z) {
  n <- nrow(X)
  m <- ncol(X)
  In <- diag(n)
  Im <- diag(m)

  z2IXXt <- z^2 * In - X %*% t(X)
  z2IXtX <- z^2 * Im - t(X) %*% X

  invz2XtX = solve(z2IXtX)
  invz2XXt = solve(z2IXXt)

  D1z = 1/n * sum(diag(z * invz2XXt))
  D2z = 1/m * sum(diag(z * invz2XtX))


  D1zp = 1/n * sum(diag(-2 * z^2 * invz2XXt^2+invz2XXt))
  D2zp = 1/m * sum(diag(-2 * z^2 * invz2XtX^2+invz2XtX))

  D1z * D2zp + D1zp * D2z

}
