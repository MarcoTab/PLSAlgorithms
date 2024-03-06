
#-----------------------------------------------------------#
# function to calculate Q of A under inner product matrix M #
#-----------------------------------------------------------#
cal_Q_fun=function(A,M){
  P_A=A%*%solve(t(A)%*%M%*%A)%*%t(A)%*%M
  if(is.null(dim(A))){
    k=length(A)
  }else{k=nrow(A)}
  return(diag(k)-P_A)
}

#-----------------------------------------------------------#
# function to calculate two blocks PLS in a single step     #
#-----------------------------------------------------------#
twoblock_pls_onestep_fun=function(mymatrix){
  # svd function sometimes returns error message
  # however, do a transpose would work
  u=NULL
  v=NULL
  
  tryCatch({svd_result = svd(mymatrix,nu=1,nv=1)
           u=svd_result$u
           v=svd_result$v},
           error = function(e) { 
             #print("this is an error")
             svd_result = svd(t(mymatrix),nv=1,nu=1)
             u<<-svd_result$v
             v<<-svd_result$u
           }
  )
  return(list(u=u,v=v))

}

#--------------------------------------#
# function to calculate two blocks pls #
#--------------------------------------#
twoblock_pls_fun=function(X,Y,d){
  i=1
  S_XY=cov(X,Y)
  S_YX=cov(Y,X)
  S_X=var(X)
  S_Y=var(Y)
  mymatrix=S_XY
  u.m=matrix(0,nrow=ncol(X),ncol=d)
  v.m=matrix(0,nrow=ncol(Y),ncol=d)
  # for(i in 1:d){
  #   twoblock_pls_onestep_result=twoblock_pls_onestep_fun(mymatrix)
  #   u.current=twoblock_pls_onestep_result$u
  #   v.current=twoblock_pls_onestep_result$v
  #   u.m[,i]=u.current
  #   v.m[,i]=v.current
  #   mymatrix=t(cal_Q_fun(u.m[,1:i],S_X))%*%S_XY%*%cal_Q_fun(v.m[,1:i],S_Y)
  # }
  twoblock_result=twoblock_pls_bymatrix_fun(SigmaX=S_X,SigmaY=S_Y,SigmaXY=S_XY,d)
  u.m=twoblock_result$u.m
  v.m=twoblock_result$v.m
  return(list(u.m=u.m,v.m=v.m))
}


#--------------------------------------#
# function to calculate two blocks pls #
#--------------------------------------#
twoblock_pls_bymatrix_fun=function(SigmaX,SigmaY,SigmaXY,d){
  i=1
  # S_XY=cov(X,Y)
  # S_YX=cov(Y,X)
  # S_X=var(X)
  # S_Y=var(Y)
  mymatrix=SigmaXY
  u.m=matrix(0,nrow=ncol(SigmaX),ncol=d)
  v.m=matrix(0,nrow=ncol(SigmaY),ncol=d)
  for(i in 1:d){
    twoblock_pls_onestep_result=twoblock_pls_onestep_fun(mymatrix)
    u.current=twoblock_pls_onestep_result$u
    v.current=twoblock_pls_onestep_result$v
    u.m[,i]=u.current
    v.m[,i]=v.current
    mymatrix=t(cal_Q_fun(u.m[,1:i],SigmaX))%*%SigmaXY%*%cal_Q_fun(v.m[,1:i],SigmaY)
  }
  return(list(u.m=u.m,v.m=v.m,singular_value.v=svd(mymatrix)$d))
}

#-------------------------------------#
# Function to calculate the dimension # 
# of two block PLS in population      #
#-------------------------------------#

cal_d_twoblock_fun=function(SigmaX,SigmaY,SigmaXY,mytol=1e-3){
  maxdim=min(nrow(SigmaX),nrow(SigmaY))
  d=0
  singular_value_sum=sum(svd(SigmaXY)$d)
  while(d<=maxdim & singular_value_sum>mytol){
    d=d+1
    twoblock_result=twoblock_pls_bymatrix_fun(SigmaX=SigmaX,SigmaY=SigmaY,SigmaXY=SigmaXY,d)
    singular_value_sum=sum(twoblock_result$singular_value.v)
  }
  return(d)
}





#----------------------------------------------------#
# function to calculate simultaneuos pls in one step #
#----------------------------------------------------#
simul_pls_onstep_fun=function(mymatrix){
  u=svd(mymatrix)$u[,1]
  return(list(u=u))
}

#----------------------------------------#
# function to calculate simultaneous pls #
#----------------------------------------#
simul_pls_fun=function(X,Y,d1,d2){
  i=1
  j=1
  S_XY=cov(X,Y)
  S_YX=cov(Y,X)
  S_X=var(X)
  S_Y=var(Y)
  mymatrix1=S_XY%*%S_YX
  mymatrix2=S_YX%*%S_XY
  u.m=matrix(0,nrow=ncol(X),ncol=d1)
  v.m=matrix(0,nrow=ncol(Y),ncol=d2)
  for(i in 1:d1){
    simul_pls_onestep_result=simul_pls_onstep_fun(mymatrix1)
    u.current=simul_pls_onestep_result$u
    u.m[,i]=u.current
    mymatrix1=t(cal_Q_fun(u.m[,1:i],S_X))%*%S_XY%*%S_YX%*%cal_Q_fun(u.m[,1:i],S_X)
  }
  if(d2==ncol(Y)){
    v.m=diag(d2)
  }else{
    for(j in 1:d2){
      simul_pls_onestep_result=simul_pls_onstep_fun(mymatrix2)
      v.current=simul_pls_onestep_result$u
      v.m[,j]=v.current
      mymatrix2=t(cal_Q_fun(v.m[,1:j],S_Y))%*%S_YX%*%S_XY%*%cal_Q_fun(v.m[,1:j],S_Y)
    }
  }
  
  
  return(list(u.m=u.m,v.m=v.m))
}


#---------------------------------------------#
# function to calculate the hat beta hat Y    #
#      and MSE for two block pls              #
#---------------------------------------------#
cal_twoblockpls_fun=function(d,X_CV_train,Y_CV_train,X_CV_test,Y_CV_test, 
                             scaleX_CV_train.v, scaleY_CV_train.v){
  X=X_CV_train
  Y=Y_CV_train
  S_XY=cov(X,Y)
  S_YX=cov(Y,X)
  S_X=var(X)
  S_Y=var(Y)
  mean_X.m=matrix(rep(colMeans(X),each=nrow(X_CV_test)),ncol=ncol(X))
  mean_Y.m=matrix(rep(colMeans(Y),each=nrow(Y_CV_test)),ncol=ncol(Y))
  
  twoblock_result=twoblock_pls_fun(X,Y,d=d)
  twoblock_hat_beta=twoblock_result$u.m%*%solve(t(twoblock_result$u.m)%*%S_X%*%twoblock_result$u.m)%*%
    t(twoblock_result$u.m)%*%S_XY%*%twoblock_result$v.m%*%t(twoblock_result$v.m)
  twoblock_hat_Y=mean_Y.m+(X_CV_test%*%(diag(1/scaleX_CV_train.v))-mean_X.m)%*%twoblock_hat_beta #### changed
  twoblock_Y_MSE=sum((Y_CV_test%*%(diag(1/scaleY_CV_train.v))-twoblock_hat_Y)^2)  ### changed
  CV_MSE_twoblock_Y_test.v<<-(Y_CV_test%*%(diag(1/scaleY_CV_train.v))-twoblock_hat_Y)^2
  cor_twoblock_Y_Yhat<<-cor(Y_CV_test%*%(diag(1/scaleY_CV_train.v)),twoblock_hat_Y)
  return(c(twoblock_Y_MSE))
}

#---------------------------------------------#
# function to calculate the hat beta hat Y    #
#      and MSE for simultaneous pls           #
#---------------------------------------------#
cal_simul_pls_fun=function(d.v,X_CV_train,Y_CV_train,X_CV_test,Y_CV_test, 
                           scaleX_CV_train.v, scaleY_CV_train.v){
  print(d.v)
  X=X_CV_train
  Y=Y_CV_train
  #Y=as.matrix(Y[,4]) ###############*****************
  S_XY=cov(X,Y)
  S_YX=cov(Y,X)
  S_X=var(X)
  S_Y=var(Y)
  d1=d.v[1]
  d2=d.v[2]
  mean_X.m=matrix(rep(colMeans(X),each=nrow(X_CV_test)),ncol=ncol(X))
  mean_Y.m=matrix(rep(colMeans(Y),each=nrow(Y_CV_test)),ncol=ncol(Y))
  
  # pls_result=plsr(Y~X,scale=F,validation="cv")
  # print(pls_result$validation)
  # pls_hat_beta=pls_result$coefficient[,,35]
  # pls_hat_Y=predict(pls_result,X_CV_test)[,,35]
  
  simul_pls_result=simul_pls_fun(X,Y,d1=d1,d2=d2)
  simul_pls_hat_beta=simul_pls_result$u.m%*%solve(t(simul_pls_result$u.m)%*%S_X%*%simul_pls_result$u.m)%*%
    t(simul_pls_result$u.m)%*%S_XY%*%simul_pls_result$v.m%*%t(simul_pls_result$v.m)
  simul_pls_hat_Y=mean_Y.m+(X_CV_test%*%(diag(1/scaleX_CV_train.v))-mean_X.m)%*%simul_pls_hat_beta
  simul_pls_Y_MSE=sum((Y_CV_test%*%(diag(1/scaleY_CV_train.v))-simul_pls_hat_Y)^2)
  CV_MSE_simul_pls_Y_test.v<<-(Y_CV_test%*%(diag(1/scaleY_CV_train.v))-simul_pls_hat_Y)^2   ###changed
  cor_simul_pls_Y_Yhat<<-cor(Y_CV_test%*%(diag(1/scaleY_CV_train.v)),simul_pls_hat_Y)    ### changed
  return(c(simul_pls_Y_MSE))
}
#------------------------------------------------------------#
#       function to use cross validation to select           #
# the dimension of two blocks PLS, simultaneous PLS and PLS1 #
#------------------------------------------------------------#
cv_fun=function(X,Y,nfold,if.Xscale,if.Yscale){
  #Y=as.matrix(Y[,4]) #####***********************
  k=nrow(X)%/%nfold # number of obs in each fold if not the last fold
  d.m=expand.grid(1:min(ncol(X),nrow(X)-k-1,20),1:min(ncol(Y),nrow(Y)-k-1,20))
  twoblock_d.v=1:min(ncol(X),min(ncol(Y),nrow(Y)-k-1),20)
  # d.m=expand.grid(1:min(ncol(X),nrow(X)-k-1),1:min(ncol(Y),nrow(Y)-k-1))
  # twoblock_d.v=1:min(ncol(X),min(ncol(Y),nrow(Y)-k-1))
  # d.m=d.m[1:3,] ######*********
  print(d.m)
  d_pls.m=expand.grid(1:min(ncol(X),nrow(X)-k-1,20),ncol(Y))
  d_pls1.m=expand.grid(1:min(ncol(X),nrow(X)-k-1,20),1)
  colnames(d.m)=c("dim_X","dim_Y")
  colnames(d_pls.m)=c("dim_X","dim_Y")
  colnames(d_pls1.m)=c("dim_X","dim_Y")
  
  #  d.m=matrix(c(35,1),ncol=2) ###########*************
  twoblock_CV_Y_MSE.m<<-matrix(0,ncol=min(ncol(X),length(twoblock_d.v)),nrow=nfold)
  simul_pls_CV_Y_MSE.m<<-matrix(0,ncol=nrow(d.m),nrow=nfold)
  pls_CV_Y_MSE.m<<-matrix(0,ncol=nrow(d_pls.m),nrow=nfold)
 # pls1_CV_Y_MSE.m<<-matrix(0,ncol=nrow(d_pls1.m),nrow=nfold)
  pls1_CV_Y_MSE.m<<-array(0,c(nfold, ncol(Y), nrow(d_pls1.m)))
  pls1_CV_onefold_result=matrix(0,nrow=ncol(Y),ncol=nrow(d_pls1.m))
  for(i in 1:nfold){
    cat(sprintf("********Below is to calculate the %d th fold********\n",i))
    CV_test_startobs=1+(i-1)*k
    if(i!=nfold){
      CV_test_endobs=k+(i-1)*k
    }else{
      CV_test_endobs=nrow(X)
    }
    X_CV_test=X[CV_test_startobs:CV_test_endobs,]
    Y_CV_test=as.matrix(Y[CV_test_startobs:CV_test_endobs,])
    X_CV_train=X[-(CV_test_startobs:CV_test_endobs),]
    Y_CV_train=as.matrix(Y[-(CV_test_startobs:CV_test_endobs),])
    if(if.Xscale){
      scaleX_CV_train.v= apply(X_CV_train,2,sd)
      X_CV_train_scaled = sweep(X_CV_train,2,scaleX_CV_train.v,'/')
    }else{
      scaleX_CV_train.v= rep(1,length=ncol(X_CV_train))
      X_CV_train_scaled = X_CV_train
    }
    if(if.Yscale){
      scaleY_CV_train.v= apply(Y_CV_train,2,sd)
      Y_CV_train_scaled = sweep(Y_CV_train,2,scaleY_CV_train.v,'/')
    }else{
      scaleY_CV_train.v= rep(1,length=ncol(Y_CV_train))
      Y_CV_train_scaled = Y_CV_train
    }
    
    twoblock_CV_onefold_result=apply(as.matrix(twoblock_d.v),1,cal_twoblockpls_fun,
                                     X_CV_train=X_CV_train_scaled,Y_CV_train=Y_CV_train_scaled,
                                     X_CV_test=X_CV_test,Y_CV_test=Y_CV_test, 
                                     scaleX_CV_train.v = scaleX_CV_train.v, 
                                     scaleY_CV_train.v = scaleY_CV_train.v)
    twoblock_CV_Y_MSE.m[i,]<<-twoblock_CV_onefold_result
    print(twoblock_CV_onefold_result)
    simul_pls_CV_onefold_result=apply(d.m,1,cal_simul_pls_fun,
                                      X_CV_train=X_CV_train_scaled,Y_CV_train=Y_CV_train_scaled,
                                      X_CV_test=X_CV_test,Y_CV_test=Y_CV_test,
                                      scaleX_CV_train.v = scaleX_CV_train.v, 
                                      scaleY_CV_train.v = scaleY_CV_train.v)
    simul_pls_CV_Y_MSE.m[i,]<<-simul_pls_CV_onefold_result
    print(simul_pls_CV_onefold_result)
    pls_CV_onefold_result=apply(d_pls.m,1,cal_simul_pls_fun,
                                X_CV_train=X_CV_train_scaled,Y_CV_train=Y_CV_train_scaled,
                                X_CV_test=X_CV_test,Y_CV_test=Y_CV_test,
                                scaleX_CV_train.v = scaleX_CV_train.v, 
                                scaleY_CV_train.v = scaleY_CV_train.v)
    pls_CV_Y_MSE.m[i,]<<-pls_CV_onefold_result
    
    for (j in 1:ncol(Y)){
#      if(j==1){temp=rep(0,nrow(d_pls1.m))}
      pls1_CV_onefold_result[j,]=apply(d_pls1.m,1,cal_simul_pls_fun,
                     X_CV_train=X_CV_train_scaled,Y_CV_train=as.matrix(Y_CV_train_scaled[,j]),
                     X_CV_test=X_CV_test,Y_CV_test=as.matrix(Y_CV_test[,j]),
                     scaleX_CV_train.v = scaleX_CV_train.v, 
                     scaleY_CV_train.v = scaleY_CV_train.v[j])

    pls1_CV_Y_MSE.m[i,j,]<<-pls1_CV_onefold_result[j,]
    }
  }
  twoblock_CV_Y_MSE.v<<-apply(twoblock_CV_Y_MSE.m,2,sum)
  d_selected=(1:min(ncol(X),ncol(Y)))[which.min(twoblock_CV_Y_MSE.v)]
  simul_pls_CV_Y_MSE.v<<-apply(simul_pls_CV_Y_MSE.m,2,sum)
  d.v_selected=d.m[which.min(simul_pls_CV_Y_MSE.v),]
  pls_CV_Y_MSE.v<<-apply(pls_CV_Y_MSE.m,2,sum)
  d_pls.v_selected=d_pls.m[which.min(pls_CV_Y_MSE.v),]
  pls1_CV_Y_MSE.v<<-apply(pls1_CV_Y_MSE.m,2:3,sum)
#  d_pls1.v_selected=d_pls1.m[which.min(pls1_CV_Y_MSE.v),]
  d_pls1.v_selected=d_pls1.m[apply(pls1_CV_Y_MSE.v, 1, which.min),]
#  d_pls1.v_selected=d_pls1.m[which.min(pls1_CV_Y_MSE.v),]
  rownames(d_pls1.v_selected)=NULL
  return(list(d=d_selected,d1=d.v_selected[1],d2=d.v_selected[2],d_pls=d_pls.v_selected,d_pls1=d_pls1.v_selected,d.m=d.m,
              twoblock_CV_Y_MSE.m=twoblock_CV_Y_MSE.m,simul_pls_CV_Y_MSE.m=simul_pls_CV_Y_MSE.m,
              pls_CV_Y_MSE.m=pls_CV_Y_MSE.m,pls1_CV_Y_MSE.m=pls1_CV_Y_MSE.m))
}


refit_entire_training_and_get_testMSE_fun=function(Xtrain,Ytrain,Xtest,Ytest,d,d1,d2,d_pls,d_pls1,
                                                   if.Xscale,if.Yscale){
  if(if.Xscale){
    scaleXtrain.v= apply(Xtrain,2,sd)
  }else{
    scaleXtrain.v= rep(1,length=ncol(Xtrain))
  }
  
  Xtrain_scaled = sweep(Xtrain,2,scaleXtrain.v,'/')
  
  if(if.Yscale){
    scaleYtrain.v= apply(Ytrain,2,sd)
  }else{
    scaleYtrain.v= rep(1,length=ncol(Ytrain))
  }
  
  Ytrain_scaled = sweep(Ytrain,2,scaleYtrain.v,'/')
  
  X=Xtrain_scaled
  Y=Ytrain_scaled
  S_XY=cov(X,Y)
  S_YX=cov(Y,X)
  S_X=var(X)
  S_Y=var(Y)
  mean_X.m=matrix(rep(colMeans(Xtrain),each=nrow(Xtest)),ncol=ncol(X))
  mean_Y.m=matrix(rep(colMeans(Ytrain),each=nrow(Ytest)),ncol=ncol(Y))
  
  twoblock_hat_beta=matrix(NA,nrow=ncol(X),ncol=ncol(Y))
  twoblock_hat_Y=Ytest*NA
  twoblock_Y_MSE.v=Ytest*NA
  tryCatch(
    {
      twoblock_result=twoblock_pls_fun(X,Y,d=d)
      twoblock_hat_beta=twoblock_result$u.m%*%solve(t(twoblock_result$u.m)%*%S_X%*%twoblock_result$u.m)%*%
          t(twoblock_result$u.m)%*%S_XY%*%twoblock_result$v.m%*%t(twoblock_result$v.m)
      twoblock_hat_Y=mean_Y.m+(Xtest-mean_X.m)%*%(diag(1/scaleXtrain.v))%*%twoblock_hat_beta%*%(diag(scaleYtrain.v)) #### changed
      twoblock_Y_MSE.v=(Ytest-twoblock_hat_Y)^2  
    },
      error=function(cond) {
      # message("Here's the original error message:")
      message(cond)
    },
      warning=function(cond) {
      # message("Here's the original warning message:")
      message(cond)
  }
  )    
  
  d1=unlist(d1)
  d2=unlist(d2)
  
  simul_pls_hat_beta=matrix(NA,nrow=ncol(X),ncol=ncol(Y))
  simul_pls_hat_Y=Ytest*NA
  simul_pls_Y_MSE.v=Ytest*NA
  
  tryCatch(
    {
      simul_pls_result=simul_pls_fun(X,Y,d1=d1,d2=d2)
      simul_pls_hat_beta=simul_pls_result$u.m%*%solve(t(simul_pls_result$u.m)%*%S_X%*%simul_pls_result$u.m)%*%
        t(simul_pls_result$u.m)%*%S_XY%*%simul_pls_result$v.m%*%t(simul_pls_result$v.m)
      simul_pls_hat_Y=mean_Y.m+(Xtest-mean_X.m)%*%(diag(1/scaleXtrain.v))%*%simul_pls_hat_beta%*%(diag(scaleYtrain.v))
      simul_pls_Y_MSE.v=(Ytest-simul_pls_hat_Y)^2
    },
  error=function(cond) {
      # message("Here's the original error message:")
      message(cond)
  },
      warning=function(cond) {
      # message("Here's the original warning message:")
      message(cond)
  }
  )    
  
  pls_hat_beta=matrix(NA,nrow=ncol(X),ncol=ncol(Y))
  pls_hat_Y=Ytest*NA
  pls_Y_MSE.v=Ytest*NA
  tryCatch(
    {
      pls_result=simul_pls_fun(X,Y,d1=unlist(d_pls[1]),d2=unlist(d_pls[2]))
      pls_hat_beta=pls_result$u.m%*%solve(t(pls_result$u.m)%*%S_X%*%pls_result$u.m)%*%
        t(pls_result$u.m)%*%S_XY%*%pls_result$v.m%*%t(pls_result$v.m)
      pls_hat_Y=mean_Y.m+(Xtest-mean_X.m)%*%(diag(1/scaleXtrain.v))%*%pls_hat_beta%*%(diag(scaleYtrain.v))
      pls_Y_MSE.v=(Ytest-pls_hat_Y)^2
    },
      error=function(cond) {
      # message("Here's the original error message:")
      message(cond)
  },
      warning=function(cond) {
      # message("Here's the original warning message:")
      message(cond)
  }
  )   
  
  
  pls1_hat_beta=matrix(NA,nrow=ncol(X),ncol=ncol(Y))
  pls1_hat_Y=Ytest*NA
  pls1_Y_MSE.v=Ytest*NA
  tryCatch(
    {
      for (j in 1:ncol(Y)){
        pls1_result=simul_pls_fun(X,as.matrix(Y[,j]),d1=unlist(d_pls1[j,1]),d2=unlist(d_pls1[j,2]))
        pls1_hat_beta[,j]=pls1_result$u.m%*%solve(t(pls1_result$u.m)%*%S_X%*%pls1_result$u.m)%*%
            t(pls1_result$u.m)%*%matrix(S_XY[,j],ncol=1)%*%pls1_result$v.m%*%t(pls1_result$v.m)
      }
    pls1_hat_Y=mean_Y.m+(Xtest-mean_X.m)%*%(diag(1/scaleXtrain.v))%*%pls1_hat_beta%*%(diag(scaleYtrain.v))
    pls1_Y_MSE.v=(Ytest-pls1_hat_Y)^2
    },
    error=function(cond) {
      # message("Here's the original error message:")
      message(cond)
    },
    warning=function(cond) {
      # message("Here's the original warning message:")
      message(cond)
    }
  )  
    return(list(twoblock_Y_MSE.v = twoblock_Y_MSE.v,
              simul_pls_Y_MSE.v = simul_pls_Y_MSE.v,
              pls_Y_MSE.v = pls_Y_MSE.v,
              pls1_Y_MSE.v = pls1_Y_MSE.v,
              twoblock_hat_Y=twoblock_hat_Y, 
              simul_pls_hat_Y=simul_pls_hat_Y,
              pls_hat_Y=pls_hat_Y,
              pls1_hat_Y=pls1_hat_Y,
              twoblock_hat_beta = twoblock_hat_beta, 
              simul_pls_hat_beta = simul_pls_hat_beta,
              pls_hat_beta = pls_hat_beta,
              pls1_hat_beta = pls1_hat_beta
    ))
}


#-----------------------------------------------#
# function to calculate leave one out MSE       #
# This is the measure in the Skagerberg paper   #
#-----------------------------------------------#

cal_LOOMSE_onefold_fun=function(i,myX,myY,d,d1,d2,d_pls,d_pls1,if.Xscale,if.Yscale){
  Xtrain_LOO=as.matrix(myX[-i,])
  Ytrain_LOO=as.matrix(myY[-i,])
  Xtest_LOO=matrix(myX[i,],nrow=1)
  Ytest_LOO=matrix(myY[i,],nrow=1)
  result_1_time_temp = refit_entire_training_and_get_testMSE_fun(Xtrain = Xtrain_LOO, Ytrain = Ytrain_LOO,
                                                                 Xtest = Xtest_LOO, Ytest = Ytest_LOO,
                                                                 d, d1, d2, d_pls, d_pls1, if.Xscale, if.Yscale)
  result_1_refit_iter_siml_pls=refit_iterative_simul_pls_fun(Xtrain = Xtrain_LOO, Ytrain = Ytrain_LOO,
                                                             Xtest = Xtest_LOO, Ytest = Ytest_LOO,d1=d1,d2=d2)
  result_1_time=list(twoblock_Y_MSE.v=(result_1_time_temp$twoblock_hat_Y-Ytest_LOO)^2,
                     simul_pls_Y_MSE.v=(result_1_time_temp$simul_pls_hat_Y-Ytest_LOO)^2,
                     iter_simul_pls_Y_MSE.v=(result_1_refit_iter_siml_pls$iter_simul_pls_hat_Y-Ytest_LOO)^2,
                     pls_Y_MSE.v=(result_1_time_temp$pls_hat_Y-Ytest_LOO)^2,
                     pls1_Y_MSE.v=(result_1_time_temp$pls1_hat_Y-Ytest_LOO)^2,
                     twoblock_hat_Y=result_1_time_temp$twoblock_hat_Y,
                     simul_pls_hat_Y=result_1_time_temp$simul_pls_hat_Y,
                     iter_simul_pls_hat_Y=result_1_refit_iter_siml_pls$iter_simul_pls_hat_Y,
                     pls_hat_Y=result_1_time_temp$pls_hat_Y,
                     pls1_hat_Y=result_1_time_temp$pls1_hat_Y)
  
  
return(result_1_time)
}


#------------------------------------------------#
# function to calculate simulatenous pls         #
# with iteration between X and Y entirely        #
#------------------------------------------------#
refit_iterative_simul_pls_fun=function(Xtrain,Ytrain,Xtest,Ytest,d1,d2
                                               ){
  # run no scaling only 
  X=Xtrain
  Y=Ytrain
  S_XY=cov(X,Y)
  S_YX=cov(Y,X)
  S_X=var(X)
  S_Y=var(Y)
  mean_X.m=matrix(rep(colMeans(Xtrain),each=nrow(Xtest)),ncol=ncol(X))
  mean_Y.m=matrix(rep(colMeans(Ytrain),each=nrow(Ytest)),ncol=ncol(Y))

  d1=unlist(d1)
  d2=unlist(d2)
  myX=X
  myY=Y
  myd1=d1
  myd2=d2
  iter_simul_pls_result=simul_pls_fun(myX,myY,d1,d2)
  myU=iter_simul_pls_result$u.m
  # iter_simul_pls_result=plsr(myY~myX,myd1,scale=F)
  # myU=iter_simul_pls_result$loadings
  myV=diag(ncol(Y))
  for(l in 1:10){
    if(l%%2 ==1){
      myU=iter_simul_pls_result$u.m
      myX=myX%*%myU
      myY=Y
      iter_simul_pls_result=simul_pls_fun(myX,myY,d1,d2)
    }else{
      myV=iter_simul_pls_result$v.m
      myY=myY%*%myV
      myX=X
      iter_simul_pls_result=simul_pls_fun(myX,myY,d1,d2)
    }
  }
  iter_simul_pls_hat_beta=myU%*%solve(t(myU)%*%S_X%*%myU)%*%
    t(myU)%*%S_XY%*%myV%*%t(myV)
  iter_simul_pls_hat_Y=mean_Y.m+(Xtest-mean_X.m)%*%iter_simul_pls_hat_beta
  iter_simul_pls_Y_MSE.v=(Ytest-iter_simul_pls_hat_Y)^2
  
  return(list(iter_simul_pls_Y_MSE.v = iter_simul_pls_Y_MSE.v,
              iter_simul_pls_hat_Y=iter_simul_pls_hat_Y,
              iter_simul_pls_hat_beta = iter_simul_pls_hat_beta))
}
