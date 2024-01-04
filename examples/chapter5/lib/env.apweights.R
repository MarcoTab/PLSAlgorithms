env.apweights <- function(X, Y, u) {

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
		
	} else if (u == r) {
	
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
		
	} else {
		
		ite <- 100
		epsilon <- 1e-5
		C <- diag(n)
		l2 <- rep(0, ite)
		Xn <- X
		Yn <- Y
		
		for (i in 1:ite) {
			m <- env(Xn, Yn, u)
			l2[i] <- m$loglik
			resid <- Yn - as.matrix(rep(1, n)) %*% t(m$mu) - Xn %*% t(m$beta)
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
		alphahat <- m$mu
		betahat <- m$beta
		Sigmahat <- m$Sigma
		Omegahat <- m$Omega
		Omega0hat <- m$Omega0
		loglik <- m$loglik
	}
 
   return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, alpha = alphahat, beta = betahat, Sigma = Sigmahat, Omega = Omegahat, Omega0 = Omega0hat, loglik = loglik,C=C1))
}
	


