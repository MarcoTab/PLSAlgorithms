
library(chemometrics)
ycal <- read.table("data/chapter4/ycal.txt", quote="\"", comment.char="")
 
ytest <- read.table("data/chapter4/ytest.txt", quote="\"", comment.char="")
 
Xcal <- read.table("data/chapter4/Xcal.txt", quote="\"", comment.char="")
 
Xtest <- read.table("data/chapter4/Xtest.txt", quote="\"", comment.char="")
 



ycal=as.vector(t(ycal))
XXcal=t(Xcal) 
XXtest=t(Xtest)

m=dim(Xcal)[1]

X1=XXcal[,seq(1,m,length=10)]

X1test=XXtest[,seq(1,m,length=10)]

X2=XXcal[,seq(1,m,length=15)]

X2test=XXtest[,seq(1,m,length=15)]

X3=XXcal[,seq(1,m,length=33)]

X3test=XXtest[,seq(1,m,length=33)]

X4=XXcal[,seq(1,m,length=50)]

X4test=XXtest[,seq(1,m,length=50)]

X5=XXcal

X5test=XXtest






gas1 <- plsr(ycal~X1, ncomp = 4,scale=TRUE)

a_10=sqrt(mean((predict(gas1, ncomp = 4, newdata = (X1test))-ytest$V1)^2))




gas2 <- plsr(ycal~X2, ncomp = 4,scale=TRUE)

a_15=sqrt(mean((predict(gas2, ncomp = 4, newdata = (X2test))-ytest$V1)^2))



gas3 <- plsr(ycal~X3, ncomp = 4,scale=TRUE)

a_33=sqrt(mean((predict(gas3, ncomp = 4, newdata = (X3test))-ytest$V1)^2))


gas4 <- plsr(ycal~X4, ncomp = 4,scale=TRUE)

a_50=sqrt(mean((predict(gas4, ncomp = 4, newdata = (X4test))-ytest$V1)^2))


gas5 <- plsr(ycal~X5, ncomp = 4,scale=TRUE)

a_101=sqrt(mean((predict(gas5, ncomp = 4, newdata = (X5test))-ytest$V1)^2))

plot(c(10,20, 33,50,101),xlab="p",c(a_10,a_15, a_33,a_50,a_101),type='b',ylab=TeX("$$\\sqrt{MSE}$$"), col="black", lwd=2, main="Figure 4.6")
