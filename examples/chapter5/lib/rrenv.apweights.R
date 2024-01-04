rrenv.apweights <- function(X, Y, u, d) {

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
	
		tmp <- env.apweights(X, Y, u)
		Gammahat <- tmp$Gamma
		Gamma0hat <- tmp$Gamma0	
		betahat <- tmp$beta
		alphahat <- tmp$alpha
		Omegahat <- tmp$Omega
		Omega0hat <- tmp$Omega0
		Sigmahat <- tmp$Sigma
		loglik <- tmp$loglik
	
	} else {
		
		ite <- 100
		epsilon <- 1e-5
		C <- diag(n)
		l2 <- rep(0, ite)
		Xn <- X
		Yn <- Y
		
		for (i in 1:ite) {
			m <- rrenv(Xn, Yn, u, d)
			l2[i] <- m$loglik
			resid <- Yn - as.matrix(rep(1, n)) %*% t(m$alpha) - Xn %*% t(m$beta)
			C1 <- diag(resid %*% solve(m$Sigma) %*% t(resid)) / r
			Xn <- diag(1 / sqrt(C1)) %*% X
			Yn <- diag(1 / sqrt(C1)) %*% Y
			if (i > 1) {
				if (abs(l2[i] - l2[i-1]) < epsilon * l2[i-1]) break			
			}
		}
		print(i)
		Gammahat <- m$Gamma
		Gamma0hat <- m$Gamma0
		alphahat <- m$alpha
		betahat <- m$beta
		Sigmahat <- m$Sigma
		Omegahat <- m$Omega
		Omega0hat <- m$Omega0
		loglik <- m$loglik
	}
 
   return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, alpha = alphahat, beta = betahat, Sigma = Sigmahat, Omega = Omegahat, Omega0 = Omega0hat, loglik = loglik,C=C1))
}
	


