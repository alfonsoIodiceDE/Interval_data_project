int_centering<-function(centers, radii, method="centers"){
library(tidyverse)
  mediac<-colMeans(C) 
  mediar<-colMeans(R) 
  a<-matrix(1,n,1) 
  mc<-a%*%mediac
  mr<-a%*%mediar
  
  if(method=="centers"){  
    
    Cov<-cov(C)
    iS<-solve(diag(diag(Cov))^0.5)
    tol=0.00001
    if(((sum(mediac)<tol)&& (abs(sum(diag(C)-p))<tol))){
      print("I dati sono standardizzati")
      Cc<-C
      Rc<-R
    }else{
      Cc<-(C-mc)%*%iS
      Rc<-(R-mr)%*%iS
    }
    vart<-(t(Cc)%*%Cc+t(Rc)%*%Rc+t(abs(Cc))%*%abs(Rc)+t(abs(Rc))%*%abs(Cc))/n
    var1<-diag(vart)
    cc=Cc%*%diag(var1^(-1/2))
    dd=Rc%*%diag(var1^(-1/2))
  }else{
    interval_min=(C-R)
    names(interval_min)=paste0(names(interval_min),"_min")
    interval_max=(C+R)
    names(interval_max)=paste0(names(interval_max),"_max")
    mean_mat=rep(1,nrow(C)) %*% t(mediac %>% as.vector()) 
    Gu=matrix(0,nrow(C),ncol(C))
    Gu[C>mean_mat]=1
    Gu[!C>mean_mat]=-1
    
    Qu = (interval_min - mean_mat)^2+
      (interval_max - mean_mat)*(interval_min - mean_mat)+(interval_max - mean_mat)^2
    
    names(Qu)=str_replace(names(Qu),pattern="_min",replacement = "")
    
    GuQu = as.matrix(Gu*sqrt(Qu))
    GuQu_t_GuQu = t(GuQu) %*% GuQu
    
    Cov= GuQu_t_GuQu/(3*nrow(C))
    iS = diag(1/sqrt(diag(Cov)))
    tol=0.00001
    if(((sum(mediac)<tol)&& (abs(sum(diag(C)-p))<tol))){
      print("I dati sono standardizzati")
      Cc<-C
      Rc<-R
    }else{
      cc<-(C-mc)%*%iS
      dd<-(R-mr)%*%iS
    }
    #Corr=Cov/sqrt(den_corr)   
    #iS = diag(1/sqrt(diag(Cov)))
  }
  
  out=list()
  out$Cov=Cov
  out$cc=cc
  out$dd=dd
  return(out)
}
