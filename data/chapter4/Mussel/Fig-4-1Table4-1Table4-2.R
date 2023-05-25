library(pls)
#library(plsr)
library(car)
# predictores L,H,S,W

musselslog <- read.csv("musselslogMinus3.txt", sep="")
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


###MARCO: esta es la figura 4.1 que es lo primero que hay que reproducir
### this is the plot of Mussel data (observed Y vs the fitted vauel for the pls fit with one
##component, we first show that leave one out choose one component)


mod <- plsr(y ~ X, ncomp = 4, scale=FALSE,validation = "LOO")
summary(mod)

#this is the picture 4.1 
 
modM <- plsr(y~X, ncomp = 1)
#pdf("yvsfitted.pdf", bg = "transparent", width = 8, height = 6)
plot(modM$fitted.values,y,main="",
     pch=19,ylab=TeX("$Y_i$"),xlab=TeX("$\\hat{Y}_i$"))
abline(0,1,col="blue",lwd=2)
dev.off()



###MARCO: esta es la tabla 4.1 Table 4.1: Coefficient estimates and 
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

printstr = "    | Estimate | S.D.\n-----------------------\n"
subs = c("₁", "₂", "₃", "₄")
for (i in 1:4) {
    betachar = paste("β", subs[i], sep="")
    printstr <- paste(printstr, betachar, "|")
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

 