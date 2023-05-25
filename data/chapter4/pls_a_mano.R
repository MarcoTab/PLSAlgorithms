#library(glasso)
#library(ridgecov)
#library(psych)


#simple is the firstr and envelope is the thirs
library(Renvlp)
estimate.beta.pls=function(X, Y, d,plsS=TRUE, plsN=TRUE,  env_PLS=TRUE, quiet=FALSE,how.much=.995)
{
  n=dim(X)[1]
  degf=n-1
  p=dim(X)[2]
  mX=apply(X,2,mean)
  if (is.null(dim(Y))){my=mean(Y)}
  else{my=apply(Y,2,mean)}
  Xc=scale(X, center=mX, scale=FALSE)
  yc=scale(Y,center=my,scale=FALSE)
  sXX=crossprod(Xc)/degf
  sXy=crossprod(Xc,yc)/degf
  syy=crossprod(yc,yc)/degf
  
  
   
  
  b.ols.pls=NULL
  b.env.pls=NULL
  if(env_PLS)
  {
    
    auxx=NULL
    
    ss=(abs(t(cov(X,Y))%*%eigen(cov(X))$vectors))
    
    if (is.null(dim(Y))){k=order(ss,decreasing=TRUE)
    for (i in 1:p){
      auxx[i]=sum(ss[k[1:i]])/sum(ss)
    }
    
    s=which.max(auxx>how.much)
    mm=k[1:s]
    }
    else{kk=ss
    st=NULL
    auxx=matrix(0,nrow=p,ncol=nrow(ss))
    for (jj in 1:nrow(ss)){
         kk[jj,]=order(ss[jj,],decreasing=TRUE)}
      for (jj in 1:nrow(ss)){
      for (i in 1:p){
        auxx[i,jj]=sum(ss[jj,kk[jj,1:i]])/sum(ss[jj,])
      }
      }
      
      st=apply((auxx>how.much),2,which.max)
      
      
      cc=NULL
      for (i in 1:nrow(kk)){
        
       cc=c(cc,c(kk[i,1:st[i]]))
      }
     mm=unique(cc)  
    }
       
     if (n>p){mm=p
    XX=Xc
    miro=xenv(XX,yc,min(d,mm))
    b.env.pls=miro$beta
    aa=lm(yc~XX-1)
    b.ols.pls=aa$coefficients
    aa=lm(yc~Xc)
    b.pc.pls=NULL
    env=miro$Gamma
    }
    
   
    if (n<p){
      XX=Xc%*%eigen(cov(X))$vectors[,mm]
      cat("eigenvectors = ", mm, "\n") 
      miro=xenv(XX,yc,min(d,mm))
      b.env.pls=t(t(miro$beta)%*%t(eigen(cov(X))$vectors[,mm]))
      aa=lm(yc~XX-1)
      b.ols.pls=t(t(aa$coefficients)%*%t(eigen(cov(X))$vectors[,mm]))
      XX=Xc%*%eigen(cov(X))$vectors[,mm[1:d]]
      aa=lm(yc~XX-1)
      b.pc.pls=t(t(aa$coefficients)%*%t(eigen(cov(X))$vectors[,mm[1:d]]))
      env=eigen(cov(X))$vectors[,mm]%*%miro$Gamma
          }
     
    
  } 
    

  
  
  b.ols=NULL
  b.env=NULL
   
    auxx=NULL
 
    ss=eigen(cov(X))$values
    if (is.null(dim(Y))){k=order(ss,decreasing=TRUE)
    for (i in 1:p){
      auxx[i]=sum(ss[k[1:i]])/sum(ss)
    }
    
    s=which.max(auxx>how.much)
    mm=k[1:s]
    }
    else{kk=ss
    st=NULL
    auxx=matrix(0,nrow=p,ncol=nrow(ss))
    for (jj in 1:nrow(ss)){
      kk[jj,]=order(ss[jj,],decreasing=TRUE)}
    for (jj in 1:nrow(ss)){
      for (i in 1:p){
        auxx[i,jj]=sum(ss[jj,kk[jj,1:i]])/sum(ss[jj,])
      }
    }
    }
    
     
     
      
     
    
    if (n>p){mm=p
    XX=Xc
    miro=xenv(XX,yc,d)
    b.env=miro$beta
    aa=lm(yc~XX-1)
    b.ols=aa$coefficients
    aa=lm(yc~Xc)
    b.pc=NULL
    
    }
    
    
    
    
    #if (n>p){XX=Xc}
    #cat("number of eigenvectors = ", s, "\n") 
    if (n<p){
      XX=Xc%*%eigen(cov(X))$vectors[,mm]
      cat("eigenvectors = ", mm, "\n") 
      miro=xenv(XX,yc,min(d,mm))
      b.env=t(t(miro$beta)%*%t(eigen(cov(X))$vectors[,mm]))
      aa=lm(yc~XX-1)
      b.ols=t(t(aa$coefficients)%*%t(eigen(cov(X))$vectors[,mm]))
      XX=Xc%*%eigen(cov(X))$vectors[,mm[1:d]]
      aa=lm(yc~XX-1)
      b.pc=t(t(aa$coefficients)%*%t(eigen(cov(X))$vectors[,mm[1:d]]))
       
    }

    
    
 
    
    
    
    
    
 
  
  
  
  
  
  
  
  if(plsS){
    # if(d==1)
    # {
    #   
    #   ## pls of order 1
    #   
    #   aux1=(sXy)%*%t(sXy)
    #   aux= eigen(aux1,symmetric=TRUE)$vectors[,1]
    #   b.plsS=aux%*%solve(t(aux)%*%sXX%*%aux)%*%t(aux)%*%(sXy)
    #   #b.pls=(eigen(sXX,symmetric=T)$values[1])^(-1)*sXy
    # }
    # if (d==2)
    # {
    #   
    #   aux1=(sXy)%*%t(sXy)
    #   v1= eigen(aux1,symmetric=TRUE)$vectors[,1]
    #   auxx=diag(p)-(as.numeric(t(v1)%*%sXX%*%v1))^(-1)*v1%*%t(v1)%*%sXX
    #   aux2=t(auxx)%*%aux1%*%auxx
    #   v2=eigen(aux2,symmetric=TRUE)$vectors[,1]
    #   H=cbind(v1,v2)
    #   b.plsS=H%*%solve(t(H)%*%sXX%*%H)%*%t(H)%*%sXy
    # }
    # 
    # if (d==3)
    # {
    #   
    #   aux1=(sXy)%*%t(sXy)
    #   v1= eigen(aux1,symmetric=TRUE)$vectors[,1]
    #   auxx=diag(p)-(as.numeric(t(v1)%*%sXX%*%v1))^(-1)*v1%*%t(v1)%*%sXX
    #   aux2=t(auxx)%*%aux1%*%auxx
    #   v2=eigen(aux2,symmetric=TRUE)$vectors[,1]
    #   
    #   H=cbind(v1,v2)
    #   auxx=diag(p)-H%*%solve(t(H)%*%sXX%*%H)%*%t(H)%*%sXX
    #   aux2=t(auxx)%*%aux1%*%auxx
    #   v3=eigen(aux2,symmetric=TRUE)$vectors[,1]
    #   
    #   H=cbind(v1,v2,v3)
    #   b.plsS=H%*%solve(t(H)%*%sXX%*%H)%*%t(H)%*%sXy
    # }
    # 
    aux1=(sXy)%*%t(sXy)
    v1= eigen(aux1,symmetric=TRUE)$vectors[,1]
    H=cbind(v1)
    for (i in 2:d){
      auxx=diag(p)-H%*%solve(t(H)%*%sXX%*%H)%*%t(H)%*%sXX
      aux2=t(auxx)%*%aux1%*%auxx
      vi=eigen(aux2,symmetric=TRUE)$vectors[,1]
      H=cbind(H,vi)
    }
    b.plsS=H%*%solve(t(H)%*%sXX%*%H)%*%t(H)%*%sXy
      
  }
  
  
   
  
  
  
  if(plsN){
    # if(d==1)
    # {
    #   
    #   ## pls of order 1
    #   
    #   aux1=(sXy)%*%solve(syy)%*%t(sXy)
    #   aux= eigen(aux1,symmetric=TRUE)$vectors[,1]
    #   b.plsN=aux%*%solve(t(aux)%*%sXX%*%aux)%*%t(aux)%*%(sXy)
    #   #b.pls=(eigen(sXX,symmetric=T)$values[1])^(-1)*sXy
    # }
    # if (d==2)
    # {
    #   
    #   aux1=(sXy)%*%solve(syy)%*%t(sXy)
    #   v1= eigen(aux1,symmetric=TRUE)$vectors[,1]
    #   auxx=diag(p)-(as.numeric(t(v1)%*%sXX%*%v1))^(-1)*v1%*%t(v1)%*%sXX
    #   aux2=t(auxx)%*%aux1%*%auxx
    #   v2=eigen(aux2,symmetric=TRUE)$vectors[,1]
    #   H=cbind(v1,v2)
    #   b.plsN=H%*%solve(t(H)%*%sXX%*%H)%*%t(H)%*%sXy
    # }
    # 
    # if (d==3)
    # {
    #   
    #   aux1=(sXy)%*%solve(syy)%*%t(sXy)
    #   v1= eigen(aux1,symmetric=TRUE)$vectors[,1]
    #   auxx=diag(p)-(as.numeric(t(v1)%*%sXX%*%v1))^(-1)*v1%*%t(v1)%*%sXX
    #   aux2=t(auxx)%*%aux1%*%auxx
    #   v2=eigen(aux2,symmetric=TRUE)$vectors[,1]
    #   
    #   H=cbind(v1,v2)
    #   auxx=diag(p)-H%*%solve(t(H)%*%sXX%*%H)%*%t(H)%*%sXX
    #   aux2=t(auxx)%*%aux1%*%auxx
    #   v3=eigen(aux2,symmetric=TRUE)$vectors[,1]
    #   
    #   H=cbind(v1,v2,v3)
    #   b.plsS=H%*%solve(t(H)%*%sXX%*%H)%*%t(H)%*%sXy
    # }
    
    aux1=(sXy)%*%solve(syy)%*%t(sXy)
    v1= eigen(aux1,symmetric=TRUE)$vectors[,1]
    H=cbind(v1)
    for (i in 2:d){
      auxx=diag(p)-H%*%solve(t(H)%*%sXX%*%H)%*%t(H)%*%sXX
      aux2=t(auxx)%*%aux1%*%auxx
      vi=eigen(aux2,symmetric=TRUE)$vectors[,1]
      H=cbind(H,vi)
    }
    b.plsN=H%*%solve(t(H)%*%sXX%*%H)%*%t(H)%*%sXy
    
  }
  
  
  
  return(cbind(b.plsS,b.plsN,b.env.pls,b.ols.pls,b.env,b.ols,b.pc))
    
}



#estimate PC and ENV in the first PLS

estimate.beta.pc=function(X, Y, d,  env_PC=TRUE, quiet=FALSE,how.much=.995)
{
  n=dim(X)[1]
  degf=n-1
  p=dim(X)[2]
  mX=apply(X,2,mean)
  if (is.null(dim(Y))){my=mean(Y)}
  else{my=apply(Y,2,mean)}
  Xc=scale(X, center=mX, scale=FALSE)
  yc=scale(Y,center=my,scale=FALSE)
  sXX=crossprod(Xc)/degf
  sXy=crossprod(Xc,yc)/degf
  syy=crossprod(yc,yc)/degf
  
  
  
  
  
  if(env_PC)
  {
    
    
    aux= eigen(cov(X),symmetric=TRUE)
    auxx=NULL
      for (i in 1:p){
          auxx[i]=sum(aux$values[1:i])
        }
      s=which.max(auxx/auxx[p]>how.much)
    
    
     
     
    
    
    
    
    if (n>p){s=p}
    cat("number of eigenvectors = ", s, "\n")   
    if (n<p){
      XX=Xc%*%eigen(cov(X))$vectors[,1:s]
    }
    cat("eigenvectors = ", 1:s, "\n") 
    dd=u.xenv(XX,yc)
    b.env.pc=t(t(xenv(XX,yc,min(d,p))$beta)%*%t(eigen(cov(X))$vectors[,1:s]))
    
    
    XXX=XX[,1:d]
    aa=lm(yc~XXX-1)
    b.ols.pc=t(t(aa$coefficients)%*%t(eigen(cov(X))$vectors[,1:d]))
    
    
  }
  
  
  
  
  return(cbind(b.env.pc,b.ols.pc))
  
}





estimate.d=function(X, Y, how.much=.995)
{
  n=dim(X)[1]
  degf=n-1
  p=dim(X)[2]
  mX=apply(X,2,mean)
  if (is.null(dim(Y))){my=mean(Y)}
  else{my=apply(Y,2,mean)}
  Xc=scale(X, center=mX, scale=FALSE)
  yc=scale(Y,center=my,scale=FALSE)
  sXX=crossprod(Xc)/degf
  sXy=crossprod(Xc,yc)/degf
  syy=crossprod(yc,yc)/degf
  
  
  
  
   
    
    # aux= eigen(cov(X),symmetric=TRUE)
    auxx=NULL
    #  for (i in 1:p){
    #     auxx[i]=sum(aux$values[1:i])
    #   }
    # s=which.max(auxx/auxx[p]>.995)
    
    
    ss=c(abs(t(cov(X,Y))%*%eigen(cov(X))$vectors))
    
    
    k=order(ss,decreasing=TRUE)
    for (i in 1:p){
      auxx[i]=sum(ss[k[1:i]])/sum(ss)
    }
    
    s=which.max(auxx>how.much)
    
    
    
    
    
    if (n>p){s=p}
    if (n>p){XX=Xc}
    cat("number of eigenvectors = ", s, "\n") 
    if (n<p){
      XX=Xc%*%eigen(cov(X))$vectors[,k[1:s]]
    }
    if (n<p){
      cat("eigenvectors = ", k[1:s], "\n") }
     dd=u.xenv(XX,yc)
     
    
     return(dd)
    
  }
