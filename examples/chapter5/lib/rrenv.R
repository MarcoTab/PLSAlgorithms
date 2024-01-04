rrenv <- function(X, Y, u, d) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
 
  
  sigY <- cov(Y) * (n - 1) / n
  sigYX <- cov(Y, X) * (n - 1) / n
  sigX <- cov(X) * (n - 1) / n
  invsigX <- chol2inv(chol(sigX))
  betaOLS <- sigYX %*% invsigX
  
  U <- tcrossprod(betaOLS, sigYX) 
  M <- sigY - U


  
	if (u == 0) {
	
		tmp <- envMU(M, U, u)
		Gammahat <- tmp$Gammahat
		Gamma0hat <- tmp$Gamma0hat
		etahat <- NULL
		Omegahat <- NULL
		Bhat <- NULL
		Omega0hat <- sigY
		alphahat <- colMeans(Y)
		betahat <- matrix(0, r, p)
		Sigmahat <- sigY
		
	} else if (u == r & d == r) {
	
		tmp <- envMU(M, U, u)
		Gammahat <- tmp$Gammahat
		Gamma0hat <- tmp$Gamma0hat
		etahat <- betaOLS
		Bhat <- diag(r)
		Omegahat <- diag(r)
		Omega0hat <- NULL
		alphahat <- colMeans(Y) - betaOLS %*% colMeans(X)
		betahat <- betaOLS
		Sigmahat <- M
		
	} else if (u == r & d < r) {
	
		tmp <- envMU(M, U, u)
		Gammahat <- tmp$Gammahat
		Gamma0hat <- tmp$Gamma0hat
		eigSx <- eigen(sigX)
		Sxnhalf <- sweep(eigSx$vectors, MARGIN = 2, 1 / sqrt(eigSx$values), '*') %*% t(eigSx$vectors) 
		eigsigY <- eigen(sigY)
		sigYnhalf <- sweep(eigsigY$vectors, MARGIN = 2, 1 / sqrt(eigsigY$values), '*') %*% t(eigsigY$vectors)
		sigYhalf <- sweep(eigsigY$vectors, MARGIN = 2, sqrt(eigsigY$values), '*') %*% t(eigsigY$vectors)
		CYX <- sigYnhalf %*% sigYX %*% Sxnhalf
		svdC <- svd(CYX)
		svdCd <- sweep(svdC$u[, 1:d], MARGIN = 2, svdC$d[1:d], '*') %*% t(svdC$v[, 1:d])
		betahat <- sigYhalf %*% svdCd %*% Sxnhalf
		alphahat <- colMeans(Y) - betahat %*% colMeans(X)
		Sigmahat <- sigYhalf %*% (diag(r) - svdCd %*% t(svdCd)) %*% sigYhalf
		Omegahat <- Sigmahat
		Omega0hat <- NULL

	} else if (d == u | d == min(p, r)) {
	
		tmp <- envMU(M, U, u)
		Gammahat <- tmp$Gammahat
		Gamma0hat <- tmp$Gamma0hat	
		etahat <- crossprod(Gammahat, betaOLS)
		betahat <- Gammahat %*% etahat
		alphahat <- colMeans(Y) - betahat %*% colMeans(X)
		Omegahat <- crossprod(Gammahat, M) %*% Gammahat
		Omega0hat <- crossprod(Gamma0hat, sigY) %*% Gamma0hat
		Sigma1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
		Sigmahat <- Sigma1 + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)	
	
	} else {
		tmp <- rrenvMU(M, U, u, d)
		Gammahat <- tmp$Gammahat
		Gamma0hat <- tmp$Gamma0hat
		eigSx <- eigen(sigX)
		Sxnhalf <- sweep(eigSx$vectors, MARGIN = 2, 1 / sqrt(eigSx$values), '*') %*% t(eigSx$vectors) 
		SgamY <- crossprod(Gammahat, sigY) %*% Gammahat
		eigSgamY <- eigen(SgamY)
		SgamYnhalf <- sweep(eigSgamY$vectors, MARGIN = 2, 1 / sqrt(eigSgamY$values), '*') %*% t(eigSgamY$vectors)
		SgamYhalf <- sweep(eigSgamY$vectors, MARGIN = 2, sqrt(eigSgamY$values), '*') %*% t(eigSgamY$vectors)
		CgamYX <- SgamYnhalf %*% t(Gammahat) %*% sigYX %*% Sxnhalf
		svdC <- svd(CgamYX)
		svdCd <- sweep(svdC$u[, 1:d], MARGIN = 2, svdC$d[1:d], '*') %*% t(svdC$v[, 1:d])
		betahat <- Gammahat %*% SgamYhalf %*% svdCd %*% Sxnhalf
		alphahat <- colMeans(Y) - betahat %*% colMeans(X)
		Omegahat <- SgamYhalf %*% (diag(u) - svdCd %*% t(svdCd)) %*% SgamYhalf
		Omega0hat <- crossprod(Gamma0hat, sigY) %*% Gamma0hat
		Sigma1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
		Sigmahat <- Sigma1 + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)

	}
	
	t1 <- eigen(Sigmahat)
	invSigmahat <- sweep(t1$vectors, MARGIN = 2, 1 / t1$values, '*') %*% t(t1$vectors)
	loglik <- -n / 2 * sum(log(t1$values)) - n / 2 * sum(diag(invSigmahat %*% (M + (betaOLS - betahat) %*% sigX%*% t(betaOLS - betahat))))
 
   return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, alpha = alphahat, beta = betahat, Sigma = Sigmahat, Omega = Omegahat, Omega0 = Omega0hat, loglik = loglik))
}
	


