rm(list=ls())
library(caret)
source("examples/chapter7/lib/lineal_with_R_more_classes.R")
source("examples/chapter7/lib/pls_matrices.R")
library(readr)
library(mltools)
library(data.table)


library(caret)

library(readr)
FTIR_Spectra_olive_oils <- read_csv("data/chapter7/FTIR_Spectra_olive_oils.csv")


Y=as.factor(t(FTIR_Spectra_olive_oils[2,2:121]))



Y_numeric=as.numeric(Y)
Y=matrix(Y_numeric,ncol=1)


xx=t(FTIR_Spectra_olive_oils[3:570,2:121])

X=matrix(0,ncol=ncol(xx),nrow=nrow(xx))
for (i in 1:ncol(xx)){
  for (j in 1:nrow(xx)){
    X[j,i]=as.numeric(xx[j,i])
  }
}


Dmax=33 

K=5


### Use this loop to choose d

error_lineal_PFC=NULL
error_lineal_pls=NULL
error_lineal_iso=NULL
set.seed(1)
for (d in 1:Dmax){

  aux=lineal.junto.small.R.more.classes(X,Y,K,d)
  
  li=as.numeric(aux$predict_lineal_pls)
  li_PFC=as.numeric(aux$predict_lineal_PFC)
  li_iso=as.numeric(aux$predict_lineal_iso)
  
  YY=as.numeric(Y)
  error_lineal_pls[d]=mean((li-YY)^2)
  error_lineal_PFC[d]=mean((li_PFC-YY)^2)
  error_lineal_iso[d]=mean((li_iso-YY)^2)
  
  
}






d_pls=which(error_lineal_pls==min(error_lineal_pls))
d_PFC=which(error_lineal_PFC==min(error_lineal_PFC))
d_iso=which(error_lineal_iso==min(error_lineal_iso))
 

aux_pls=lineal.junto.small.R.more.classes(X,Y,nrow(X),d_pls[1])$predict_lineal_pls
aux_PFC=lineal.junto.small.R.more.classes(X,Y,nrow(X),d_PFC[1])$predict_lineal_PFC
aux_iso=lineal.junto.small.R.more.classes(X,Y,nrow(X),d_iso[1])$predict_lineal_iso



# For table 7.1, the second row, first 3 columns

col1 = confusionMatrix(as.factor(aux_pls),as.factor(Y))$overall[1]
col2 = confusionMatrix(as.factor(aux_PFC),as.factor(Y ))$overall[1]
col3 = confusionMatrix(as.factor(aux_iso),as.factor(Y ))$overall[1]

cat("Table 7.1\nDataset   | PLS  | ISO | PFC\nOlive Oil |", round(col1,2), "|", round(col2,2), "|", round(col3,2), "\n")

### Figure 7.4

#######pls 
AUX=pls2_nipals(X ,Y ,d)$P[,1:d_pls]

datos_proy=X%*%AUX


grayPalette <- c("#555555", "#AAAAAA", "#717171", "#8D8D8D")


dfp_proj <- as.data.frame(datos_proy)
origins <- c("Greece", "Italy", "Portugal", "Spain")
dfp_proj$origins <- origins[Y] 

ggplot(dfp_proj, aes(x = V1, y = V2, shape = origins, color = origins)) + 
  geom_point(alpha = 0.7, size = 4) +
  theme_bw() +
  scale_color_manual(values=grayPalette) +
  labs(x = "First PLS projection",
       y = "Second PLS projection") +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 




dfp_proj <- as.data.frame(datos_proy)
origins <- c("Greece", "Italy", "Portugal", "Spain")
dfp_proj$origins <- origins[Y] 

ggplot(dfp_proj, aes(x = V1, y = V2, shape = origins, color = origins)) + 
  geom_point(alpha = 0.7, size = 4) +
  theme_bw() +
  labs(x = "First PLS projection",
       y = "Second PLS projection") +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 





###### iso

reg=lm(X ~Y )
vec=eigen(cov(reg$residuals))$vectors[,1:d_iso]
datos_proy_iso=X%*%vec


dfp_proj_iso <- as.data.frame(datos_proy_iso)
origins <- c("Greece", "Italy", "Portugal", "Spain")
dfp_proj_iso$origins <- origins[Y] 

ggplot(dfp_proj_iso, aes(x = V1, y = V2, shape = origins, color = origins)) + 
  geom_point(alpha = 0.7, size = 4) +
  theme_bw() +
  labs(x = "First ISO projection",
       y = "Second ISO projection") +
  scale_color_manual(values=grayPalette) +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 





dfp_proj_iso <- as.data.frame(datos_proy_iso)
origins <- c("Greece", "Italy", "Portugal", "Spain")
dfp_proj_iso$origins <- origins[Y] 

ggplot(dfp_proj_iso, aes(x = V1, y = V2, shape = origins, color = origins)) + 
  geom_point(alpha = 0.7, size = 4) +
  theme_bw() +
  labs(x = "First ISO projection",
       y = "Second ISO projection") +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 




#### pfc



AUX=pls2_nipals(X ,Y ,d)$P[,1:d_PFC]

datos_proy=X%*%AUX







s0 <- dr(Y~datos_proy)
datos_proy_pfc=datos_proy%*%s0$evectors[,1:2] 




dfp_proj_pfc <- as.data.frame(datos_proy_pfc)
origins <- c("Greece", "Italy", "Portugal", "Spain")
dfp_proj_pfc$origins <- origins[Y] 

ggplot(dfp_proj_pfc, aes(x = Dir1, y = Dir2, shape = origins, color = origins)) + 
  geom_point(alpha = 0.7, size = 4) +
  theme_bw() +
  labs(x = "First PFC projection",
       y = "Second PFC projection") +
  scale_color_manual(values=grayPalette) +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 



dfp_proj_pfc <- as.data.frame(datos_proy_pfc)
origins <- c("Greece", "Italy", "Portugal", "Spain")
dfp_proj_pfc$origins <- origins[Y] 

ggplot(dfp_proj_pfc, aes(x = Dir1, y = Dir2, shape = origins, color = origins)) + 
  geom_point(alpha = 0.7, size = 4) +
  theme_bw() +
  labs(x = "First PFC projection",
       y = "Second PFC projection") +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size=16),
        legend.text = element_text(size = 13)
  ) 






 