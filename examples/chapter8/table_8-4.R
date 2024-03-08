### tomo los datos...
rm(list=ls())
library(caret)
  
library(readr)

 
FTIR_Spectra_instant_coffee <- read.csv("data/chapter8/coffee.csv", header=FALSE)
 
source("examples/chapter8/lib/quadratic_with_cross_validation.R")
source("examples/chapter8/lib/pls_matrices.R")
Y=as.numeric(t(FTIR_Spectra_instant_coffee[2,2:57]))
 
Y=matrix(Y,ncol=1)
xx=t(FTIR_Spectra_instant_coffee[4:289,2:57])

 
 

X=matrix(0,ncol=ncol(xx),nrow=nrow(xx))
for (i in 1:ncol(xx)){
  for (j in 1:nrow(xx)){
    X[j,i]=as.numeric(xx[j,i])
  }
}



Dmax=12 

K=10
h=seq(0,1,.2) 

quad_1=matrix(0,ncol=length(h),nrow=Dmax)
quad_2=matrix(0,ncol=length(h),nrow=Dmax)

set.seed(1)

cat("Fitting coffee data\n")
for (d in 1:Dmax){
  for (s in 1:length(h)){
    cat(d, ", ", s, " / ", Dmax,", ", length(h), '\n', sep="")
    
    alp=h[s]
    aux=quadratic.discrimination(X,Y,K,d,alp)
    li_1=aux$predict_quad_1
    li_2=aux$predict_quad_2
    
    quad_1[d,s]=mean((li_1-Y)^2)
    quad_2[d,s]=mean((li_2-Y)^2)
  }
}


M1=which(quad_1 == min(quad_1), arr.ind = TRUE)
M2=which(quad_2 == min(quad_2), arr.ind = TRUE)
#minimums
#quad 1 d=4, alp=0.8
#quad_2 d=2, alp =  1 


####Table 8.4: Olive oil and Coffee data: 
###Estimates of the correct classification rates (%) from leave one out cross validation. Results in columns 4–8 are as described in Table 7.1.

d=M1[1]
alp=M1[2]

aux=quadratic.discrimination(X,Y,nrow(X),d,alp)
cofanq1 <- confusionMatrix(as.factor(aux$predict_quad_1),as.factor(Y))$overall[1]


d=M2[1]
alp=M2[2]

aux=quadratic.discrimination(X,Y,nrow(X),d,alp)
cofanq2 = confusionMatrix(as.factor(aux$predict_quad_2),as.factor(Y))$overall[1]





oils <- read_csv("data/chapter8/FTIR_Spectra_olive_oils.csv", col_names = FALSE)

FTIR_Spectra_instant_coffee=oils


xx=t(FTIR_Spectra_instant_coffee[4:573,2:121])

X=matrix(0,ncol=ncol(xx),nrow=nrow(xx))
for (i in 1:ncol(xx)){
  for (j in 1:nrow(xx)){
    X[j,i]=as.numeric(xx[j,i])
  }
}



Y=as.factor(t(FTIR_Spectra_instant_coffee[3,2:121]))
Y=as.numeric(Y)
Y=matrix(Y,ncol=1)

Dmax=14
K=nrow(X)

h=seq(0,1,.2) 

quad_1=matrix(0,ncol=length(h),nrow=Dmax)
quad_2=matrix(0,ncol=length(h),nrow=Dmax)

cat("Fitting olive oil data\n")
for (d in 1:Dmax){
  for (s in 1:length(h)){
    cat(d, ", ", s, " / ", Dmax,", ", length(h), '\n', sep="")
    
    alp=h[s]
    
    aux=quadratic.discrimination(X,Y,K,d,alp)
    
    li_1=aux$predict_quad_1
    li_2=aux$predict_quad_2
    
    
    quad_1[d,s]=mean((li_1-Y)^2)
    quad_2[d,s]=mean((li_2-Y)^2)
    
    
  }
  
}



M1=which(quad_1 == min(quad_1), arr.ind = TRUE)
M2=which(quad_2 == min(quad_2), arr.ind = TRUE)

#quad_2 d= 10, 1 
#quad 1 d= 13 1

aux=quadratic.discrimination(X,Y,nrow(X),M1[1],h[M1[2]])

library(matrixcalc)

oanq1 = confusionMatrix(as.factor(aux$predict_quad_1),as.factor(Y))$overall[1]

aux=quadratic.discrimination(X,Y,nrow(X),M2[1],h[M2[2]])

oanq2<- confusionMatrix(as.factor(aux$predict_quad_2),as.factor(Y))$overall[1]


cat("Table 8.4\n")
cat("Dataset   | AN-Q1 | AN-Q2\n")
cat("----------|-------|------\n")
cat("Coffee    |", 100*cofanq1,"|", 100*cofanq2,"\n")
cat("Olive Oil |", 100*oanq1, "|", 100*oanq2, "\n")
