rm(list=ls())
library(caret)
library(readr)
library(matrixcalc)
library(chemometrics)
FTIR_Spectra_instant_coffee <-  read_csv("data/chapter7/coffee.csv")
 
source("examples/chapter7/pls_matrices.R")
source("examples/chapter7/lineal_with_R.R")
Y=(t(FTIR_Spectra_instant_coffee[2,2:57]))
Y=as.factor(Y)

Y=as.numeric(Y)
Y=Y-1

Y=matrix(Y,ncol=1)
xx=t(FTIR_Spectra_instant_coffee[3:288,2:57])

X=matrix(0,ncol=ncol(xx),nrow=nrow(xx))
for (i in 1:ncol(xx)){
  for (j in 1:nrow(xx)){
    X[j,i]=as.numeric(xx[j,i])
  }
}


### Choose d with k=10 folds
Dmax=15 

K=10

   
linear_error_pfc=NULL
linear_error_pls=NULL
linear_error_iso=NULL
set.seed(11)
for (d in 1:Dmax){
    
    Y=c(Y)
    aux=lineal.junto.small.R(X,Y,K,d)
    #aux=lineal.junto.small.R.more.classes(X,Y,K,d)
     
    li=(aux$predict_lineal_pls)
    li_pfc=(aux$predict_lineal_PFC)
    li_iso=(aux$predict_lineal_iso)
    
    linear_error_pls[d]=mean((li-Y)^2)
    linear_error_pfc[d]=mean((li_pfc-Y)^2)
    linear_error_iso[d]=mean((li_iso-Y)^2)
    
    
  }
 



d_pls=which(linear_error_pls==min(linear_error_pls))
d_pfc=which(linear_error_pfc==min(linear_error_pfc))
d_iso=which(linear_error_iso==min(linear_error_iso))

###leave one out error for table 7.1

# Calculate the LOO error with these selected d's for the three methods.
aux_pls=lineal.junto.small.R(X,Y,nrow(X),d_pls[1])$predict_lineal_pls
aux_PFC=lineal.junto.small.R(X,Y,nrow(X),d_pfc[1])$predict_lineal_PFC
aux_iso=lineal.junto.small.R(X,Y,nrow(X),d_iso[1])$predict_lineal_iso

## LOESS Error (actual numbers for Table 7.1 for coffee)
plserr <- as.numeric(confusionMatrix(as.factor(aux_pls-1),as.factor(Y))$overall[1])
isoerr <- as.numeric(confusionMatrix(as.factor(aux_iso-1),as.factor(Y ))$overall[1])
pls.pfcerr <- as.numeric(confusionMatrix(as.factor(aux_PFC-1),as.factor(Y ))$overall[1])


cat("Dataset | PLS | ISO | PLS+PFC\n")
cat("Coffee  |", round(100.00*plserr,2), "|", round(100.00*isoerr,2), "|", round(100.00*pls.pfcerr,2))


#### Figure 7.2
Dmax=15
K=nrow(X)
AA=NULL
for (d in 1:Dmax){
  
  aux_pls= lineal.junto.small.R(X,Y,nrow(X),d)$predict_lineal_pls
 AA[d]=confusionMatrix(as.factor(aux_pls-1),as.factor(Y))$overall[1]
}

###### Figure 7.2

XX=1:Dmax
grayPalette <- c("#555555", "#AAAAAA", "#717171", "#8D8D8D")
A<- data.frame(x=XX,y=AA)
 
ggplot2::ggplot(A, aes(x = x, y = y)) + 
  geom_point(alpha = 0.8, size = 4) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values=grayPalette) +
  labs(x = "Number of components",
       y = "Accuracy") +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=16)
  ) 


library(mosaic)
library(ggplot2)

# PFC
d=3
AUX=pls2_nipals(X ,Y ,d)$P[,1:d]

datos_proy=X%*%AUX

## Figure 7.3a
dfp_proj <- as.data.frame(datos_proy)
dfp_proj$variety <- ifelse(dfp$Y==0, "robusta", "arabica")

ggplot(dfp_proj, aes(x = V1, y = V2, shape = variety, color = variety)) + 
  geom_point(alpha = 0.8, size = 4) +
  theme_bw() +
  labs(x = "First PLS projection",
       y = "Second PLS projection") +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  )

# Figure 7.3 b
s0 <- dr(Y~datos_proy)
datos_proy_pfc=datos_proy%*%s0$evectors[,1] 

dfp <- data.frame(X = datos_proy_pfc, Y = Y)
dfp$variety <- ifelse(dfp$Y==0, "robusta", "arabica")
dfp$alpha <- ifelse(dfp$Y == 0, 0.5, 1)
ggplot(dfp, aes(x = X, fill = variety, alpha = alpha)) + 
  geom_histogram(color = "black") +
  theme_bw() +
  labs(x = "PFC projection",
       y = "Number of cases") +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  )  + 
  scale_alpha(guide = 'none') 


#figure 7.3.c
#iso
d=Dmax
reg=lm(X ~Y )
vec=eigen(cov(reg$residuals))$vectors[,1:d]
datos_proy_iso=X%*%vec

dfp_proj_iso <- as.data.frame(datos_proy_iso)
dfp_proj_iso$variety <- ifelse(dfp$Y==0, "robusta", "arabica")

ggplot(dfp_proj_iso, aes(x = V1, y = V2, shape = variety, color = variety)) + 
  geom_point(alpha = 0.8, size = 4) +
  theme_bw() +
  labs(x = "First ISO projection",
       y = "Second ISO projection") +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 