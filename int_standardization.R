int_standardization<-function(x=NULL,centers, radii, std_meth="centers"){
  n <- dim(centers)[1] 
  p <- dim(centers)[2] 
  
  library(tidyverse)
  if(is.null(x)&&is.null(centers)){
    print("data must be provided either in intervals, or in centers and radii")
  }
  
  if(is.null(x)){
    
    C <- centers
    R <- radii
  }else{
    n <- nrow(x)
    pp <- ncol(x)
    x_min <- x[,seq(from=1,to=pp,by=2)]
    x_max <- x[,-seq(from=1,to=pp,by=2)]
    C <- (x_min+x_max)/2
    R <- (x_max-x_min)/2
  }
  
  mediac <- colMeans(C) 
  mediar <- colMeans(R) 
  a <- matrix(1,n,1) 
  mc <- a%*%mediac
  mr <- a%*%mediar
  
  if(std_meth=="centers"){  
    
    Cov <- cov(C)
    iS <- solve(diag(diag(Cov))^0.5)
    tol <- 0.00001
    if(((sum(mediac)<tol)&& (abs(sum(diag(C)-p))<tol))){
      print("I dati sono standardizzati")
      Cc <- C
      Rc <- R
    }else{
      Cc <- (C-mc)%*%iS
      Rc <- (R-mr)%*%iS
    }
    
    glob_var <- (t(Cc)%*%Cc+t(Rc)%*%Rc+(abs(t(Cc)%*%Rc))+(abs(t(Rc)%*%Cc)))/n
    var1 <- diag(glob_var)
    cc <- Cc%*%diag(var1^(-1/2))
    dd <- Rc%*%diag(var1^(-1/2))
    
    den_corr <- diag(Cov)%*%t(diag(Cov))
    Corr <- Cov/sqrt(den_corr)
  }else{
    interval_min <- (C-R)
    names(interval_min) <- paste0(names(interval_min),"_min")
    interval_max <- (C+R)
    names(interval_max) <- paste0(names(interval_max),"_max")
    mean_mat <- rep(1,nrow(C)) %*% t(mediac %>% as.vector()) 
    Gu <- matrix(0,nrow(C),ncol(C))
    Gu[C>mean_mat] <- 1
    Gu[!C>mean_mat] <- -1
    
    Qu <- (interval_min - mean_mat)^2+
      (interval_max - mean_mat)*(interval_min - mean_mat)+(interval_max - mean_mat)^2
    
    names(Qu) <- str_replace(names(Qu),pattern="_min",replacement = "")
    
    GuQu <- as.matrix(Gu*sqrt(Qu))
    GuQu_t_GuQu <- t(GuQu) %*% GuQu
    
    Cov <- GuQu_t_GuQu/(3*nrow(C))
    iS <- diag(1/sqrt(diag(Cov)))
    tol <- 0.00001
    if(((sum(mediac)<tol)&& (abs(sum(diag(C)-p))<tol))){
      print("I dati sono standardizzati")
      cc <- C
      dd <- R
    }else{
      cc <- (C-mc)%*%iS
      dd <- (R-mr)%*%iS
      den_corr <- diag(Cov)%*%t(diag(Cov))
      Corr <- Cov/sqrt(den_corr)
      glob_var <- (t(cc)%*%cc + t(dd)%*%dd + abs(t(cc)%*%dd) + abs(t(dd)%*%cc))/n
    }
  }
  
  out <- list()
  out$Cov <- Cov
  out$cc <- cc
  out$dd <- dd
  out$Corr <- Corr
  out$glob_var <- glob_var
  return(out)
}
