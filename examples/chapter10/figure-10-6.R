#
# Figure 10.6
#
# Rönkkö, M., Mcintosh, C., Antonakis, J., & Edwards, J. R. (2016). Partial least squares path 
# modeling: Time for some serious second thoughts. Journal of Operations Management, forthcoming.
# DOI: 10.1016/j.jom.2016.05.002
rm(list=ls())
library(Renvlp)
library(pls)
library(CCA)
library(matrixpls)
library(lavaan)
library(parallel)
library(pracma)
library(MASS)
library(psych)
 
source("rrenv.R")
source("rrr.R")
source("GE.R")
source("envMU.R")
####functions

library(ks)



Sqrt.Matrix<-function(M){
  if (!is.matrix(M)) stop("The argument is not a matrix")
  if (!identical(dim(M)[1], dim(M)[2])) stop("The matrix is not square")
  #if ( !(det(M) > 0))stop("The matrix is not nonnegative #definite")
  if ((dim(M)[1]==1) & (dim(M)[2]==1)) out<- as.matrix(sqrt(M))
  else {svd.M<-svd(M); Vec.M<-svd.M$u; Val.M<-diag(c(svd.M$d)); out<-Vec.M %*% sqrt(Val.M) %*% t(Vec.M)}
  return(out)}

#######################################      Figure 2     #######################################

# Set up parallel processing for the simulation

options(mc.cores = detectCores())

#set.seed(1)
#it is necessary to run with 100, 300, 1000
N <- 100
p=9
Lraw1=c(.8,-.7,60,rep(1,p),rep(-1,p))
Lraw2=c(1,.3,-0.59/60,rep(-1,p),rep(1,p))
Lraw=cbind(Lraw1,Lraw2)
a1=5
b1=.1

 
Lambda <- cbind(rbind(Lraw,matrix(rep(0,2*dim(Lraw)[1]),ncol=2)),rbind(matrix(rep(0,2*dim(Lraw)[1]),ncol=2),Lraw))
Lambda0<-(diag(rep(1,4*p+6))-Lambda%*%solve(t(Lambda)%*%Lambda)%*%t(Lambda))

auxLraw=.1*(Lraw1+9*Lraw2)



MODEL <- "
X =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + x21
Y =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10 + y11 + y12 + y13 + y14 + y15 + y16 + y17 + y18 + y19 + y20 + y21
Y ~ X
"



# Vary the population beta

BETAS <- seq(200/1500, 1000/1500,20/1500)

nrep=10
nose=array(NA,c(length(BETAS),5,nrep))


set.seed(30)
M=cbind(rbind(auxLraw,rep(0,21)),rbind(rep(0,21),auxLraw))

for (h in 1:nrep){
  print("h")
  
  for (j in 1:length(BETAS)){
    Beta=BETAS[j]
    
    eta <- mvrnorm(N,c(0,0),matrix(c(1,Beta,Beta,1),2,2))
    # data <- eta %*% t(Lambda) + rnorm(6*N)+ mvrnorm(n = 100, rep(0,6), diag(rep(1,6))-.1*Lambda%*%t(Lambda))
    data <- eta %*%(M)   +  mvrnorm(N, rep(0,4*p+6), a1*Lambda%*%solve(t(Lambda)%*%Lambda)%*%t(Lambda)+b1*Lambda0)
    #   data <- eta %*% t(Lambda) 
    colnames(data) <- c("x1","x2","x3","x4",
                        "x5","x6","x7","x8","x9","x10",
                        "x11","x12","x13","x14","x15",
                        "x16","x17","x18","x19","x20",
                        "x21","y1","y2","y3","y4",
                        "y5","y6","y7","y8","y9","y10",
                        "y11","y12","y13","y14","y15",
                        "y16","y17","y18","y19","y20",
                        "y21")
    
    
    # ML SEM estimates
    # ML SEM estimates
    # CB|SEM ML SEM estimates
    ml <- sem(MODEL,data)
    
     
    # MATRIXPLS Regression with PLS Mode A composits
    pls_them1 <- matrixpls(cor(data), MODEL)
    
    # IDEAL Correlation between ideal composites (equivalent to standardized regression)
    ideal <- cor(lm(eta[,1] ~ data[,1:(2*p+3)])$fitted,
                 lm(eta[,2] ~ data[,(2*p+4):(4*p+6)])$fitted)*N/(N-1)
    
    
    ###PLS TRUE  
    X=data[,1:(2*p+3)]
    Y=data[,(2*p+4):(4*p+6)]
    
    
    X= scale(X,center = TRUE, scale = FALSE)
    Y= scale(Y,center = TRUE, scale = FALSE)
    
    
    A=X%*%plsr(Y~X,ncomp=2)$projection[,1:2]
    B=Y%*%plsr(X~Y,ncomp=2)$projection[,1:2]
    
    
    
    A=A%*%sqrtm(cov(A))$Binv
    B=B%*%sqrtm(cov(B))$Binv
    
    au1=svd(cov(A,B))
    au2=svd(cov(B,A))
    
    pls_true=abs(cor(A%*%au1$u[,1],B%*%au2$u[,1]))
    
    
    
    
     
       
    
    
    
    #ENV
    
    
    X=data[,1:(2*p+3)]
    Y=data[,(2*p+4):(4*p+6)]
    
    
    X= scale(X,center = TRUE, scale = FALSE)
    Y= scale(Y,center = TRUE, scale = FALSE)
    
    
    GammaY=env(X, Y, 2, asy = TRUE)$Gamma
    
    
    GammaX=xenv(X, Y, 2, asy = TRUE)$Gamma
    
    GammaX=Lambda[1:(2*p+3),1:2]
    GammaY=Lambda[1:(2*p+3),1:2]
    nreps=10
    for (k in 1:nreps){
      Xnew=X%*%GammaX
      Ynew=Y%*%GammaY
      GammaX=xenv(X, Ynew, 2, asy = TRUE)$Gamma
      
      GammaY=env(Xnew, Y, 2, asy = TRUE)$Gamma
    }
    
    
    Ynew=Y%*%GammaY
    Xnew=X%*%GammaX 
    
    
    
    
    
    
    
    BetaYdadoX=env(X, Ynew, 2, asy = TRUE)$beta
    BetaXdadoY=env(Xnew, Y, 2, asy = TRUE)$beta
    
    au1=svd(cov(Xnew,Ynew))
    au2=svd(cov(Ynew,Xnew))
    
    
    res=(sum(diag(BetaYdadoX%*%BetaXdadoY)))^(1/2)
    
    res=abs(cor(Xnew%*%au1$u[,1],Ynew%*%au2$u[,1]))
    
     
    
    
    nose[j,,h]=c(pls_them1[1], ideal,inspect(ml,"std")$beta[2],res,
                   pls_true )
  }
}




estimates=apply(nose,1:2,function(x){mean(x,na.rm=TRUE)})


######new version of 

 


 


############finals
 
ideal_1000=estimates[,2]
MATRIXPLS_1000=estimates[,1]
ENVELOPE_EB_SEM_1000= estimates[,4]
 
CBS_SEM_1000=abs(estimates[,3])
 
EB_PLS_SEM_1000=estimates[,5]



 
ideal_300=estimates[,2]
MATRIXPLS_300=estimates[,1]
ENVELOPE_EB_SEM_300= estimates[,4]
 
CBS_SEM_300=abs(estimates[,3])
 
EB_PLS_SEM_300=estimates[,5]


5
ideal_100=estimates[,2]
MATRIXPLS_100=estimates[,1]
ENVELOPE_EB_SEM_100= estimates[,4]
 
CBS_SEM_100=abs(estimates[,3])
 
EB_PLS_SEM_100=estimates[,5]


library(latex2exp)
###ultimos
df <- data.frame(
  x1 = ideal_100, y1 = MATRIXPLS_100, y2 = MATRIXPLS_300,
  y3 = MATRIXPLS_1000, y4 = ENVELOPE_EB_SEM_100,y5=ENVELOPE_EB_SEM_300,
  y6=ENVELOPE_EB_SEM_1000,y7=EB_PLS_SEM_100,y8=EB_PLS_SEM_300, y9=EB_PLS_SEM_1000,
  x2=ideal_300,x3=ideal_1000,y10=CBS_SEM_100,y11=CBS_SEM_300,y12=CBS_SEM_1000,
  B=BETAS
)

 

###########Figure 10.6

pdf(file="FIGURE-10-6.pdf",height=8,width=10)

old_par <- par(mfrow = c(4,3), oma = c(4, 8, 3, 1), 
               mar = c(4, 2, 1, 1), las = 1 )
library(latex2exp)
plot(df$x1, df$y1  ,
     xlab = "", ylab = TeX("$\\widehat{\\Psi}$"))

abline(c(0,1))
plot(df$x2, df$y2 ,
     xlab = "", ylab = "")
abline(c(0,1))
plot(df$x3, df$y3  ,
     xlab = "", ylab = "")
abline(c(0,1))

plot(df$x1, df$y4,
     xlab = "", ylab = TeX("$\\widehat{\\Psi}$"))
abline(c(0,1))
plot(df$x2, df$y5,
     xlab = "", ylab = "")
abline(c(0,1))
plot(df$x3, df$y6,
     xlab = "", ylab = "")
abline(c(0,1))

plot(df$x1, df$y7,
     xlab = " ", ylab = " ")
title(ylab=TeX("${\\Psi}$"),xlab=TeX("${\\Psi}$"), line=3,cex.lab=1.5)
abline(c(0,1))
plot(df$x2, df$y8,
     xlab =  "", ylab = "")
title(xlab=TeX("${\\Psi}$"), line=3,cex.lab=1.5)
abline(c(0,1))
plot(df$x3, df$y9,
     xlab = "", ylab = "")
title(xlab=TeX("${\\Psi}$"), line=3,cex.lab=1.5)
abline(c(0,1))

plot(df$B, df$y10  ,
     xlab = "",ylab="")
title(ylab=TeX("$\\widehat{cor}( xi, eta))$"),xlab=TeX("${cor}( xi, eta)$"), line=3,cex.lab=1.5)
abline(c(0,1))
plot(df$B, df$y11 ,
     xlab = "", ylab = "")
title( xlab=TeX("${cor}( xi, eta)$"), line=3,cex.lab=1.5)
abline(c(0,1))
plot(df$B, df$y12  ,
     xlab = "", ylab = "")
title( xlab=TeX("${cor}( xi, eta)$"), line=3,cex.lab=1.5)
abline(c(0,1))




mtext(TeX("$\\widehat{cor}( xi, eta)$"), side=2, cex=1, 
      outer = TRUE, at=0.15, line = 0.6,las=0)





 



mtext(TeX("$\\widehat{\\Psi}$"), side=2, cex=1, 
      outer = TRUE, at=0.9, line = 0.6,las=0)


mtext(TeX("$\\widehat{\\Psi}$"), side=2, cex=1, at=0.65,
      outer = TRUE,   line = 0.6,las=0)

mtext(TeX("$\\widehat{\\Psi}$"), side=2, cex=1, 
      outer = TRUE, at=0.4, line = 0.6,las=0)

 

mtext("N = 100", side=3, cex=1.1, outer = TRUE, at = 0.17)
mtext("N = 300", side=3, cex=1.1, outer = TRUE, at = 0.50)
mtext("N = 1000", side=3, cex=1.1, outer = TRUE, at = 0.84)




mtext("Matrixpls", side=2, cex=0.8 , outer = TRUE, at = 0.9, adj = 1.7)
mtext("ENV", side=2, cex=0.8 , outer = TRUE, at = 0.65, adj = 3.3)
mtext("PLS", side=2, cex=0.8 , outer = TRUE, at = 0.4, adj = 3.5)
mtext("CB | SEM", side=2, cex=0.8 , outer = TRUE, at =0.15 , adj = 1.6)

par(old_par)
dev.off()
