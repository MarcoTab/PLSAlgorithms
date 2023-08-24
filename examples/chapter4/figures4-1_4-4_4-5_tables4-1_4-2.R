library(pls)
library(car)
# predictors L,H,S,W

musselslog <- read.csv("data/chapter4/musselslogMinus3.txt", sep="")
 attach(musselslog)
library(latex2exp)


y <- logM
Xp <- musselslog[,-3]
X <- as.matrix(Xp)
p <- dim(X)[2]
n <- length(y)
m <- n-1

y1 <- y-mean(y)
X1 <- scale(X,center=TRUE,scale=FALSE)


### this is the plot of Mussel data (observed Y vs the fitted vauel for the pls fit with one
##component, we first show that leave one out choose one component)


mod <- plsr(y ~ X, ncomp = 4, scale=FALSE,validation = "LOO")
summary(mod)

#this is the Figure 4.1 
 
modM <- plsr(y~X, ncomp = 1)
#pdf("yvsfitted.pdf", bg = "transparent", width = 8, height = 6)
plot(modM$fitted.values,y,main="",
     pch=19,ylab=TeX("$Y_i$"),xlab=TeX("$\\hat{Y}_i$"))
abline(0,1,col="blue",lwd=2)
dev.off()



###Table 4.1: Coefficient estimates and 
###corresponding asymptotic standard deviations (S.D.) 
###from three fits of the mussels’ data: Ordinary least squares, (OLS), 
###partial least squares with one component (PLS), and envelope with one 
###dimension (ENV). Results for OLS and ENV are from Cook (2018).

 
modM <- plsr(y~X, ncomp = 1)
 

## estimates
betapls <- modM$coefficients[,,1]
 
## compute SD using the formulas from 
## Envelopes and partial least squares regression
## R. D. Cook, I. S. Helland, Z. Su
## Pages 851-877
## Number of pages 27
## Journal of the Royal Statistical Society. Series B: Statistical Methodology
## Volume 75
## Issue number5
## Date Published - Nov 2013


shat <- cov(X,y)
sts <- as.numeric((t(shat)%*%shat))
sigma.0hat <- qr.Q(qr(shat),complete=TRUE)[,2:p] # px(p-1)
SIhat <- cov(X)
Omegahat <- as.numeric(t(shat)%*%SIhat%*%shat/sts)
Omega0hat <- t(sigma.0hat)%*%SIhat%*%sigma.0hat # (p-1)x(p-1)
SIhat2 <- (Omegahat/sts)*(shat%*%t(shat)) + sigma.0hat%*%Omega0hat%*%t(sigma.0hat)
sigmaY2hat <- var(y)

sds <- c()

#for beta1
XN=c(1,0,0,0)
Vhat <- as.numeric((Omegahat^(-2)/m)*t(XN)%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%XN)
sds <- append(sds, round(sqrt(Vhat), 4))

#for beta2
XN=c(0,1,0,0)
Vhat <- as.numeric((Omegahat^(-2)/m)*t(XN)%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%XN)
sds <- append(sds, round(sqrt(Vhat), 4))

#for beta3
XN=c(0,0,1,0)
Vhat <- as.numeric((Omegahat^(-2)/m)*t(XN)%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%XN)
sds <- append(sds, round(sqrt(Vhat), 4))

#for beta4
XN=c(0,0,0,1)
Vhat <- as.numeric((Omegahat^(-2)/m)*t(XN)%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%XN)
sds <- append(sds, round(sqrt(Vhat), 4))

betapls <- round(betapls,3)

# Print beta tables
printstr = "   | Estimate | S.D.\n-----------------------\n"

for (i in 1:4) {
    printstr <- paste(printstr, "β |")
    printstr <- paste(printstr, betapls[i], "   |")
    printstr <- paste(printstr, sds[i], "\n")
}
print("Table 4.1: PLS")
cat(printstr)

#####Table 4.2: Mussels’ muscles: Estimates of the covariance matrix 
#####ΣX from the envelope, standard and PLS (q = 1) fits. 

SIhat2 <- (Omegahat/sts)*(shat%*%t(shat)) + sigma.0hat%*%Omega0hat%*%t(sigma.0hat)

print("Table 4.2: PLS")
print(round(SIhat2,3))





###################################################################
#Figure 4.4
###################################################################

### WARNING
# What follows is a very time intensive simulation. Estimate ~4-6 hours

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
#############################
#############################
#simulation to get Figure 4.4
#############################
#############################
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









 