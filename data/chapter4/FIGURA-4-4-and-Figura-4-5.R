###################################################################
#This is the script to reproduce Figure 1 from the paper
#Asymptotic distribution of one-component partial least
#squares regression estimators in high dimension
#by Jeronimo Basa, Dennis Cook, Liliana Forzani and Miguel Marcos
###################################################################


library(psych)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(MASS)
library(glasso)



######Necessary function 

estimate.beta = function(X, y, pls=TRUE, quiet=FALSE){
  
  
  ###########################
  #Compute the beta coefficient for the regression 
  #y= E(y) + beta^T (X-E(X)) + error
  #y in R, X in R^p, Sigma_x=Cov(X), sigma_xy=cov(X,y)
  #assuming that the PLS model with d=1 is true, i.e.
  #beta = sigma_xy (sigma_xy^T Sigma_x sigma_xy)^{-1} (sigma_xy^T sigma_xy)
  ###########################
  
  
  n <- dim(X)[1]
  degf <- n
  p <- dim(X)[2]
  mX <- apply(X,2,mean) # media por columna
  my <- mean(y)
  Xc <- scale(X, center=mX, scale=FALSE)
  yc <- y-my
  sXX <- crossprod(Xc)/degf
  sXy <- crossprod(Xc,yc)/degf
  syy <- sum(yc^2)/degf
  if(pls){
    aux1 <- t(sXy)%*%sXy
    aux2 <- t(sXy)%*%sXX%*%sXy
    b.pls <- as.numeric(aux2)^(-1)*as.numeric(aux1)*sXy
  }
  return(b.pls)
}
###########################
###########################
#simulation to get Figure 1
###########################
###########################
set.seed(180)
reps <- 500
p.vec <- 2^(3:11)
pm <- p.vec[length(p.vec)]
pow <- 1/2
s <- pow
a <- 1
ovec <- c(1, rep(0, pm-1)) 
if(pow > 0){
  ovec=c(rep(1,  round(p.vec[1]^pow)  ), rep(0, p.vec[1]-round(p.vec[1]^pow)) )
  for(k in 2:length(p.vec)){
    b <- length(ovec)
    p <- p.vec[k]
    cnz <- sum(ovec != 0)
    tnz <- round(p^pow)
    nnz <- tnz-cnz
    newvec <- c(rep(1, nnz), rep(0, b-nnz) )
    ovec <- c(ovec, newvec)
  }
}
tau <- 1/2
pe.pls1 <- matrix(NA, nrow=length(p.vec), ncol=reps)
pe.pls2 <- matrix(NA, nrow=length(p.vec), ncol=reps)
nm <- pm/2; n.vec <- p.vec/2

for(i in 1:reps){
  cat("rep = ", i, "\n")
  sigma.xym <- ovec*rnorm(pm) 
  sigma.0 <- qr.Q(qr(sigma.xym),complete=TRUE)[,2:pm]
  for(k in 1:length(p.vec)){
    p <- p.vec[k]
    n <- p/2
    m <- n-1
    Cp1 <- p
    Cp2 <- p
    sigma.xymp <- sigma.xym[1:p]
    sTs <- as.numeric(crossprod(sigma.xymp,sigma.xymp))
    Omega <- sTs
    sigmaY2 <- tau^2 + sTs*Omega^(-1)
    sigma.0 <- qr.Q(qr(sigma.xymp),complete=TRUE)[,2:p]
    Omega0 <- diag(p-1)
    Sigma <- (Omega/sTs) * tcrossprod(sigma.xymp,sigma.xymp) + sigma.0%*%tcrossprod(Omega0,sigma.0)
    out.eigen <- eigen(Sigma, symmetric=TRUE)
    Sigma.sqrt <- out.eigen$vec%*%diag(out.eigen$val^0.5)%*%t(out.eigen$vec)
    X <- matrix(rnorm(n*p), nrow=n, ncol=p)%*%Sigma.sqrt
    true.beta <- Omega^(-1)*sigma.xymp
    y <- rnorm(n,X%*%true.beta, tau)
    Xnewc <- sigma.xymp
    b.pls <- estimate.beta(X=X, y=y, quiet=FALSE)
    bias <- (tau^2-sTs/Omega)*(Cp1/(m*sTs)) - sigmaY2*Cp1^2/(m^2*Omega*sTs) - sigmaY2*Cp2/(m*Omega*sTs)
    V <- as.numeric((Omega^(-2)/n)*t(Xnewc)%*%(Sigma*sigmaY2-sigma.xymp%*%t(sigma.xymp))%*%Xnewc)
    pe.pls1[k,i] <- V^(-1/2)*(t(b.pls) - t(true.beta))%*%Xnewc
    pe.pls2[k,i] <- V^(-1/2)*(t(b.pls) - t(true.beta)*(1+bias))%*%Xnewc
  }
}

#####FIGURE 4.4
nc <- 10
lim <- c(-7,7)

pdf("Hist.pdf", bg = "transparent")
par(mfrow=c(3,3))

hist(pe.pls1[1,], col="green", xlab="", 
     ylab="",nclass = nc, xlim = lim,ylim=c(0,0.4),freq=FALSE,main="n = 4, p = 8")
# Second with add=T to plot on top
hist(pe.pls2[1,pe.pls2[1,]<6], xlim=c(0,300), col=rgb(0,0,1,0.5), add=T, nclass = nc,freq=FALSE)
curve(dnorm(x, mean=0,sd=1), add=T, col="blue")

hist(pe.pls1[2,], col="green", xlab="", 
     ylab="",nclass = nc, xlim = lim,ylim=c(0,0.45),freq=FALSE,main="n = 8, p = 16")
# Second with add=T to plot on top
hist(pe.pls2[2,pe.pls2[2,]<6], xlim=c(0,300), col=rgb(0,0,1,0.5), add=T, nclass = nc,freq=FALSE)
curve(dnorm(x, mean=0,sd=1), add=T, col="blue")

hist(pe.pls1[3,], col="green", xlab="", 
     ylab="",nclass = nc, xlim = lim,ylim=c(0,0.4),freq=FALSE,main="n = 16, p = 32")
# Second with add=T to plot on top
hist(pe.pls2[3,pe.pls2[3,]<6], xlim=c(0,300), col=rgb(0,0,1,0.5), add=T, nclass = nc,freq=FALSE)
curve(dnorm(x, mean=0,sd=1), add=T, col="blue")

hist(pe.pls1[4,], col="green", xlab="", 
     ylab="",nclass = nc, xlim = lim,ylim=c(0,0.4),freq=FALSE,main="n = 32, p = 64")
# Second with add=T to plot on top
hist(pe.pls2[4,pe.pls2[4,]<6], xlim=c(0,300), col=rgb(0,0,1,0.5), add=T, nclass = nc,freq=FALSE)
curve(dnorm(x, mean=0,sd=1), add=T, col="blue")

hist(pe.pls1[5,], col="green", xlab="", 
     ylab="",nclass = nc, xlim = lim,ylim=c(0,0.4),freq=FALSE,main="n = 64, p = 128")
# Second with add=T to plot on top
hist(pe.pls2[5,], xlim=c(0,300), col=rgb(0,0,1,0.5), add=T, nclass = nc,freq=FALSE)
curve(dnorm(x, mean=0,sd=1), add=T, col="blue")

hist(pe.pls1[6,], col="green", xlab="", 
     ylab="",nclass = nc, xlim = lim,ylim=c(0,0.4),freq=FALSE,main="n = 128, p = 256")
# Second with add=T to plot on top
hist(pe.pls2[6,], xlim=c(0,300), col=rgb(0,0,1,0.5), add=T, nclass = nc,freq=FALSE)
curve(dnorm(x, mean=0,sd=1), add=T, col="blue")

hist(pe.pls1[7,], col="green", xlab="", 
     ylab="",nclass = nc, xlim = lim,ylim=c(0,0.4),freq=FALSE,main="n = 256, p = 512")
# Second with add=T to plot on top
hist(pe.pls2[7,], xlim=c(0,300), col=rgb(0,0,1,0.5), add=T, nclass = nc,freq=FALSE)
curve(dnorm(x, mean=0,sd=1), add=T, col="blue")

hist(pe.pls1[8,], col="green", xlab="", 
     ylab="",nclass = nc, xlim = lim,ylim=c(0,0.4),freq=FALSE,main="n = 512, p = 1024")
# Second with add=T to plot on top
hist(pe.pls2[8,], xlim=c(0,300), col=rgb(0,0,1,0.5), add=T, nclass = nc,freq=FALSE)
curve(dnorm(x, mean=0,sd=1), add=T, col="blue")

hist(pe.pls1[9,], col="green", xlab="", 
     ylab="",nclass = nc, xlim = lim,ylim=c(0,0.4),freq=FALSE,main="n = 1024, p = 2048")
# Second with add=T to plot on top
hist(pe.pls2[9,], xlim=c(0,300), col=rgb(0,0,1,0.5), add=T, nclass = nc,freq=FALSE)
curve(dnorm(x, mean=0,sd=1), add=T, col="blue")
dev.off()



## FIGURE 4-5


A=matrix(0,ncol=2,nrow=9)
B=matrix(0,ncol=2,nrow=9)

par(mfrow=c(1,2))
for (j in 1:9){
  A[j,]=c(median(pe.pls1[j,]),sd(pe.pls1[j,]))
  B[j,]=c(median(pe.pls2[j,]),sd(pe.pls2[j,pe.pls2[j,]<15]))
}

plot(log(p.vec),xlab="log p",A[,1],type='b',ylim=c(-3,10),ylab="Mean", col="red", lwd=2)
lines(log(p.vec),B[,1],col="blue", lwd=2, type='b')
abline(c(0,0))

plot(log(p.vec),xlab="log p",A[,2],type='b',ylim=c(0.5,4.25), ylab="Standard deviation", col="red", lwd=2)
lines(log(p.vec),B[,2], col="blue", lwd=2, type='b')
abline(c(1,0))

