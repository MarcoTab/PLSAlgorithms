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
pdf("yvsfitted.pdf", bg = "transparent", width = 8, height = 6)
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
 
## compute SD using the formular from 
## Envelopes and partial least squares regression
## R. D. Cook, I. S. Helland, Z. Su
## Pages (from-to)851-877
## Number of pages27
## JournalJournal of the Royal Statistical Society. Series B: Statistical Methodology
## Volume75
## Issue number5
## StatePublished - Nov 2013


shat <- cov(X,y)
sts <- as.numeric((t(shat)%*%shat))
sigma.0hat <- qr.Q(qr(shat),complete=TRUE)[,2:p] # px(p-1)
SIhat <- cov(X)
Omegahat <- as.numeric(t(shat)%*%SIhat%*%shat/sts)
Omega0hat <- t(sigma.0hat)%*%SIhat%*%sigma.0hat # (p-1)x(p-1)
SIhat2 <- (Omegahat/sts)*(shat%*%t(shat)) + sigma.0hat%*%Omega0hat%*%t(sigma.0hat)

#for beta1
XN=c(1,0,0,0)
Vhat <- as.numeric((Omegahat^(-2)/m)*t(XN)%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%XN)

#for beta2
XN=c(0,1,0,0)
Vhat <- as.numeric((Omegahat^(-2)/m)*t(XN)%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%XN)

#for beta3
XN=c(0,0,1,0)
Vhat <- as.numeric((Omegahat^(-2)/m)*t(XN)%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%XN)

#for beta4
XN=c(0,0,0,1)
Vhat <- as.numeric((Omegahat^(-2)/m)*t(XN)%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%XN)



 

predicciones <- NULL
IC1.vec <- matrix(data = NA, nrow = n, ncol = 2)
IC2.vec <- matrix(data = NA, nrow = n, ncol = 2)

for (i in 1:n){
  XN <- X[i,]
  XN <- as.numeric(XN)
  X1 <- X[-i,]
  YN <- y[i]
  Y1 <- y[-i]
  centroy <- mean(Y1)
  centrox <- apply(X1,2,mean)
  X1 <- scale(X1, center = TRUE, scale=FALSE)
  Y1 <- Y1 - centroy
  XN <- XN - centrox
  mod <- plsr(Y1~X1, 1)
  betahat <- mod$coefficients[,1,1]
  YhatN <- as.numeric(betahat%*%XN) + centroy
  shat <- cov(X1,Y1)
  sts <- as.numeric((t(shat)%*%shat))
  sigma.0hat <- qr.Q(qr(shat),complete=TRUE)[,2:p] # px(p-1)
  SIhat <- cov(X1)
  Omegahat <- as.numeric(t(shat)%*%SIhat%*%shat/sts)
  Omega0hat <- t(sigma.0hat)%*%SIhat%*%sigma.0hat # (p-1)x(p-1)
  SIhat2 <- (Omegahat/sts)*(shat%*%t(shat)) + sigma.0hat%*%Omega0hat%*%t(sigma.0hat)
  sigmaY2hat <- var(Y1)
  #Vhat1 <- as.numeric((Omegahat^(-2)/m)*t(XN)%*%(SIhat*sigmaY2hat - shat%*%t(shat))%*%XN)
  Vhat <- as.numeric((Omegahat^(-2)/m)*t(XN)%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%XN)
  residuo <- Y1 - as.numeric(X1%*%betahat)
  tauhat2 <- var(residuo)
  predicciones[i] <- YhatN
  #IC1 <- YhatN + c(-1,1)*Vhat^(1/2)*qnorm(0.975)
  #IC1.vec[i,] <- IC1
  IC2 <- YhatN + c(-1,1)*(Vhat+tauhat2)^(1/2)*qnorm(0.975)
  IC2.vec[i,] <- IC2
  cat("n = ", i, "\n")
}




plot(1:n, predicciones, xlab = "index", ylim = c(1,5.5), pch = 18, ylab = "")
#axis(1, at = seq(1, n, by = 10), las=1, labels = TRUE)
for (i in 1:n) {
  segments(i,IC2.vec[i,1],i,IC2.vec[i,2])
}
legend("topleft", 
       legend = c("Observed value", "Prediction value", "Prediction interval"), 
       col = c("red", "black", "black"), 
       pch = c(18,18,NA),
       lty = c(NA, NA, 1),
       lwd = c(NA,NA,2),
       bty = "n", 
       pt.cex = 1.5, 
       cex = 1, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.01, 0.01))
points(1:n, y, col = 'red', pch = 18)


# PLS vs OLS
alpha <- 0.05
modM <- plsr(y~X, ncomp = 1)
modO <- lm(y~X)
summary(modO)
betapls <- modM$coefficients[,,1]
betaols <- modO$coefficients[2:5]

CI <- confint(modO)
plot(1:p, modO$coefficients[2:(p+1)], xlab = "", pch = 18, ylab = "",ylim = c(-1,2),xaxt='n')
axis(1, at = seq(1, 4, by = 1), las=1, labels = TRUE)
for (i in 1:p) {
  segments(i,CI[(i+1),1],i,CI[(i+1),2])
}

o1 <- c(1,0,0,0)
o2 <- c(0,1,0,0)
o3 <- c(0,0,1,0)
o4 <- c(0,0,0,1)

CIbeta1 <- CIPLS1(y,X,o1,alpha=0.05)
CIbeta2 <- CIPLS1(y,X,o2,alpha=0.05)
CIbeta3 <- CIPLS1(y,X,o3,alpha=0.05)
CIbeta4 <- CIPLS1(y,X,o4,alpha=0.05)

betaevp1 <- 0.141
betaevp2 <- 0.154
betaevp3 <- 0.625
betaevp4 <- 0.206
betaevp <- c(betaevp1,betaevp2,betaevp3,betaevp4)

CIevp1 <- betaevp1 + c(-1,1)*qnorm(1-alpha/2)*0.0052
CIevp2 <- betaevp2 + c(-1,1)*qnorm(1-alpha/2)*0.0056
CIevp3 <- betaevp3 + c(-1,1)*qnorm(1-alpha/2)*0.0194
CIevp4 <- betaevp4 + c(-1,1)*qnorm(1-alpha/2)*0.0073



sep <- (1:p)+0.01
sep2 <- (1:p)+0.03
pdf("CIbetas.pdf", height=14, width=18)
plot(1:p, betaols, xlab = "", pch = 20, ylab = "",ylim=c(-1,2),xaxt="n")
par(new=TRUE)
points(sep,betapls,pch=20,col="lightgreen")
points(sep2,betaevp,pch=20,col="pink")
#plot(sep, betapls, xlab = "", pch = 20, ylab = "")
#axis(1, at = seq(1, 4, by = 1), las=1, labels = TRUE)
for (i in 1:p) {
  segments(i,CI[(i+1),1],i,CI[(i+1),2],col="red")
}
segments(sep[1],CIbeta1[1],sep[1],CIbeta1[2],col="darkblue",lwd=1)
segments(sep[2],CIbeta2[1],sep[2],CIbeta2[2],col="darkblue",lwd=1)
segments(sep[3],CIbeta3[1],sep[3],CIbeta3[2],col="darkblue",lwd=1)
segments(sep[4],CIbeta4[1],sep[4],CIbeta4[2],col="darkblue",lwd=1)

segments(sep2[1],CIevp1[1],sep2[1],CIevp1[2],col="darkgreen",lwd=1)
segments(sep2[2],CIevp2[1],sep2[2],CIevp2[2],col="darkgreen",lwd=1)
segments(sep2[3],CIevp3[1],sep2[3],CIevp3[2],col="darkgreen",lwd=1)
segments(sep2[4],CIevp4[1],sep2[4],CIevp4[2],col="darkgreen",lwd=1)


legend("topright", 
       legend = c("OLS estimation", "PLS estimation", "Envelope estimation", "OLS Confidence interval","PLS Confidence interval","Envelope Confidence interval"), 
       col = c("black","lightgreen","pink","red","darkblue","darkgreen"), 
       pch = c(20,20,20,NA,NA,NA),
       lty = c(NA,NA,NA,1,1,1),
       lwd = c(NA,NA,NA,2,2,2),
       bty = "n", 
       pt.cex = 3, 
       cex = 1.5, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.0001, 0.0001))
dev.off()

s <- 0
for (k in 1:length(y)){
  if (IC2.vec[k,1] <= y[k] & y[k]<=IC2.vec[k,2]){s <- s+1}
  else{s <- s}  
}
s/length(y)*100


pdf("plsvsols.pdf", bg = "transparent", width = 8, height = 6)
scatterplot(modM$fitted.values, modO$fitted.values,boxplots = FALSE,
            smooth = FALSE,pch=20,ylab="",xlab="",
            main="Valores predichos por PLS vs Valores predichos por OLS", col = "black",grid=FALSE)
abline(0,1,col="blue",lwd=2)
dev.off()

