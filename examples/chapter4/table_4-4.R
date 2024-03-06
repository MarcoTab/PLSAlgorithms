rm(list=ls())
library(pls)
library(xtable)
library(MASS)
library(psych)
library(latex2exp)

set.seed(10)
###############


###############################
#This is the script to reproduce Table 2 (table 4.4 in the book) from the paper
#Asymptotic distribution of one-component partial least
#squares regression estimators in high dimension
#by Jeronimo Basa, Dennis Cook, Liliana Forzani and Miguel Marcos
###############################


# First simulation
# Results for the third and fourth columns

#Generate data ###############
# 

p.vec <- 2^(8:10)
pmax <- tail(p.vec, n=1)
reps <- 1000
n <- floor(pmax^(1/2))
pow <- 3/4
s <- pow
ovec <- c(1, rep(0, pmax-1)) 

if(pow > 0){
  ovec=c(rep(1,round(p.vec[1]^pow)), rep(0, p.vec[1]-round(p.vec[1]^pow)))
  for(k in 2:length(p.vec)){
    b <- length(ovec)
    p <- p.vec[k]
    cnz <- sum(ovec != 0)
    tnz <- round(p^pow)
    nnz <- tnz-cnz
    newvec <- c(rep(1, nnz), rep(0, b-nnz))
    ovec <- c(ovec, newvec)
  }
}
sigmaxy <- ovec*rnorm(p)
sigma.0 <- qr.Q(qr(sigmaxy),complete=TRUE)[,2:p]
sTs <- as.numeric(t(sigmaxy)%*%sigmaxy)
Omega <- sTs
tau <- 1
sigmaY2 <- tau^2 + sTs*Omega^(-1)
Omega0 <- diag(p-1)
Sigma <- (sigmaxy%*%t(sigmaxy)/sTs)*Omega + sigma.0%*%Omega0%*%t(sigma.0)

out.eigen <- eigen(Sigma, symmetric=TRUE)
Sigma.sqrt <- out.eigen$vec%*%diag(out.eigen$val^0.5)%*%t(out.eigen$vec)
Xbig <- array(rep(NA, n*p*reps), dim=c(n, p, reps))
Ybig <- array(rep(NA, n*reps), dim=c(n, reps))
true.beta <- Omega^(-1)*sigmaxy
for (j in 1:reps){
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)%*%Sigma.sqrt
  y <- rnorm(n,X%*%true.beta, tau)
  Xbig[,,j] <- X
  Ybig[,j] <- y
}

###############
nmax <- dim(Xbig)[1]
pmax <- dim(Xbig)[2]
k <- dim(Xbig)[3]
pvec <- 2^(8:10)
nvec <- floor(pvec^(1/2))
XN <- mvrnorm(1, mu = rep(0,pmax), Sigma = diag(pmax))
cover1 <- matrix(NaN,nrow=length(pvec),ncol=k)
cover2 <- matrix(NaN,nrow=length(pvec),ncol=k)
alfa <- 0.05
z <- qnorm(1-alfa/2)
###############

###############
for (i in 1:length(pvec)){
  n <- nvec[i]
  p <- pvec[i]
  cat("p = ", p, "n = ", n, "\n")
  m <- n-1
  G <- XN[1:p]
  Omega <- as.numeric(t(sigmaxy[1:p])%*%sigmaxy[1:p])
  beta <- Omega^(-1)*sigmaxy[1:p]
  true.val <- as.numeric(t(beta)%*%G)
  Cp1 <- Cp2 <- p-1
  sts <- as.numeric(t(sigmaxy[1:p])%*%sigmaxy[1:p])
  bias <- -((sts/Omega - tau^2)*(Cp1/(m*sts)) + sigmaY2*Cp1^2/(m^2*Omega*sts) + sigmaY2*Cp2/(m*Omega*sts))
  for (j in 1:k){
    X <- Xbig[1:n, 1:p, j]
    y <- Ybig[1:n, j]
    mX <- apply(X,2,mean)
    X <- scale(X, center = mX, scale=FALSE)
    y <- y - mean(y)
    mod <- plsr(y~X, 1)
    betahat <- as.vector(mod$coefficients[,,1])
    pred <- as.numeric(crossprod(betahat,G))
    shat <- cov(X,y)
    SIhat <- cov(X)
    stshat <- as.numeric(crossprod(shat,shat))
    
    sigmaY2hat <- var(y)
    sigma.0hat <- qr.Q(qr(shat),complete=TRUE)[,2:p]
    Omegahat <- as.numeric(t(shat)%*%SIhat%*%shat/stshat)
    Omega0hat <- t(sigma.0hat)%*%SIhat%*%sigma.0hat
    
    SIhat2 <- (Omegahat/stshat)*(shat%*%t(shat)) + sigma.0hat%*%Omega0hat%*%t(sigma.0hat)
    Vhat <- as.numeric((Omegahat^(-2)/n)*t(G)%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%G)
    IC1 <- pred + c(-1,1)*z*Vhat^(1/2)
    
    cover1[i,j] <- (IC1[1] <= true.val) & (true.val <= IC1[2])
    cover2[i,j] <- (IC1[1] <= true.val*(1+bias)) & (true.val*(1+bias) <= IC1[2])
  }
  coverage1 <- sprintf("%.3f",apply(cover1,1,mean)[i])
  coverage2 <- sprintf("%.3f",apply(cover2,1,mean)[i])
  print(paste("n =", n, "p =", p, "3rd:", coverage1))
  print(paste("n =", n, "p =", p, "4th:", coverage2))
}
###################################################
###PRINTING THE RESULTS FOR COLUMNS THIRD AND FORTH
###################################################
library(gridExtra)
library(grid)
A=matrix(cbind(nvec,pvec,apply(cover1,1,mean),apply(cover2,1,mean)),nrow=3)

df <- data.frame(A)

colnames(df) <- c('n','p',"beta^T G" ,'(1+b) beta^T G')

tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
grid.table(df, theme=tt)
#############################  

# Second simulation
# Results for the fifth, sixth and seventh columns
# Uncomment the next two lines if you wish to run only the second simulation consistently
rm(list=ls())
set.seed(10)
#Generate data ###############
# 
p.vec <- 2^(8:10)
pmax <- tail(p.vec, n=1)
reps <- 1000
n <- floor(pmax^(1/2))
pow <- 3/4
s <- pow
ovec <- c(1, rep(0, pmax-1)) 

if(pow > 0){
  ovec=c(rep(1,round(p.vec[1]^pow)), rep(0, p.vec[1]-round(p.vec[1]^pow)))
  for(k in 2:length(p.vec)){
    b <- length(ovec)
    p <- p.vec[k]
    cnz <- sum(ovec != 0)
    tnz <- round(p^pow)
    nnz <- tnz-cnz
    newvec <- c(rep(1, nnz), rep(0, b-nnz))
    ovec <- c(ovec, newvec)
  }
}

sigmaxy <- ovec*rnorm(p)
sigma.0 <- qr.Q(qr(sigmaxy),complete=TRUE)[,2:p]
sTs <- as.numeric(t(sigmaxy)%*%sigmaxy)
Omega <- sTs
tau <- 1/2
sigmaY2 <- tau^2 + sTs*Omega^(-1)
Omega0 <- diag(p-1)
Sigma <- (sigmaxy%*%t(sigmaxy)/sTs)*Omega + sigma.0%*%Omega0%*%t(sigma.0)
out.eigen <- eigen(Sigma, symmetric=TRUE)
Sigma.sqrt <- out.eigen$vec%*%diag(out.eigen$val^0.5)%*%t(out.eigen$vec)
Xbig <- array(rep(NA, n*p*reps), dim=c(n, p, reps))
Ybig <- array(rep(NA, n*reps), dim=c(n, reps))
true.beta <- Omega^(-1)*sigmaxy
for (j in 1:reps){
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)%*%Sigma.sqrt
  y <- rnorm(n,X%*%true.beta, tau)
  Xbig[,,j] <- X
  Ybig[,j] <- y
}

###############
nmax <- dim(Xbig)[1]
pmax <- dim(Xbig)[2]
k <- dim(Xbig)[3]
pvec <- 2^(8:10)
nvec <- floor(pvec^(1/2))
cover1 <- NULL
cover2 <- NULL
cover3 <- NULL
alfa <- 0.05
z <- qnorm(1-alfa/2)
###############

#write(t(Xbig[,,1]), "myXbig.csv", sep=",", ncolumns=pmax)
#write(t(Ybig), "myYbig.csv", sep=',', ncolumns=1)
#write(t(sigmaxy), "mysigmaxy.csv", sep=',', ncolumns=1)

cover1 <- matrix(NaN,nrow=length(pvec),ncol=k)
cover2 <- matrix(NaN,nrow=length(pvec),ncol=k)
cover3 <- matrix(NaN,nrow=length(pvec),ncol=k)
###############
for (i in 1:length(pvec)){
  n <- nvec[i]
  p <- pvec[i]
  cat("p = ", p, "n = ", n, '\n')
  m <- n-1
  Omega <- as.numeric(t(sigmaxy[1:p])%*%sigmaxy[1:p])
  G <- sigmaxy[1:p]
  beta <- Omega^(-1)*sigmaxy[1:p]
  true.val <- as.numeric(t(beta)%*%G)
  Cp1 <- Cp2 <- p-1
  sts <- as.numeric(t(sigmaxy[1:p])%*%sigmaxy[1:p])
  bias <- -((sts/Omega - tau^2)*(Cp1/(m*sts)) + sigmaY2*Cp1^2/(m^2*Omega*sts) + sigmaY2*Cp2/(m*Omega*sts))
  for (j in 1:k){
    X <- Xbig[1:n, 1:p, j]
    y <- Ybig[1:n, j]
    mX <- apply(X,2,mean)
    X <- scale(X, center = mX, scale=FALSE)
    y <- y - mean(y)
    mod <- plsr(y~X, 1)
    betahat <- as.vector(mod$coefficients[,,1])
    pred <- as.numeric(crossprod(betahat,G))
    
    shat <- cov(X,y)
    SIhat <- cov(X)
    stshat <- as.numeric(crossprod(shat,shat))
    
    sigmaY2hat <- var(y)
    sigma.0hat <- qr.Q(qr(shat),complete=TRUE)[,2:p]
    Omegahat <- as.numeric(t(shat)%*%SIhat%*%shat/stshat)
    Omega0hat <- t(sigma.0hat)%*%SIhat%*%sigma.0hat
    Cp1hat <- tr(Omega0hat)
    Cp2hat <- tr(Omega0hat%*%Omega0hat)
    tau2hat <- var(mod$residuals)
    
    biashat <- -((stshat/Omegahat - tau2hat)*(Cp1hat/(m*stshat)) + sigmaY2hat*Cp1hat^2/(m^2*Omegahat*stshat) + sigmaY2hat*Cp2hat/(m*Omegahat*stshat))
    V <- as.numeric((Omega^(-2)/n)*t(G)%*%(Sigma[1:p,1:p]*sigmaY2 - sigmaxy[1:p]%*%t(sigmaxy[1:p]))%*%G)
    IC1 <- pred + c(-1,1)*z*V^(1/2)
    IC2 <- IC1/(1+biashat)
    #print(V)
    #print(IC1)
    #print(IC2)
    
    cover1[i,j] <- (IC1[1] <= true.val) & (true.val <= IC1[2])
    cover2[i,j] <- (IC1[1] <= true.val*(1+bias)) & (true.val*(1+bias) <= IC1[2])
    cover3[i,j] <- (IC2[1] <= true.val) & (true.val <= IC2[2])
  }
  coverage1 <- sprintf("%.3f",apply(cover1,1,mean)[i])
  coverage2 <- sprintf("%.3f",apply(cover2,1,mean)[i])
  coverage3 <- sprintf("%.3f",apply(cover3,1,mean)[i])
  print(paste("n =", n, "p =", p, "5th:", coverage1))
  print(paste("n =", n, "p =", p, "6th:", coverage2))
  print(paste("n =", n, "p =", p, "7th:", coverage3))
}

#############################

#####################################################################
###PRINTING THE RESULTS FOR COLUMNS FIFTH SIXTH AND SEVENTH
####################################################################
library(gridExtra)
library(grid)
A=matrix(cbind(nvec,pvec,apply(cover1,1,mean),apply(cover2,1,mean),apply(cover3,1,mean)),nrow=3)

df <- data.frame(A)

colnames(df) <- c('n','p',"beta^T G" ,'(1+b) beta^T G')

tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
grid.table(df, theme=tt)
############################# 