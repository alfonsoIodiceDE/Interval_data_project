interval_PCA<- function(x=NULL,centers=NULL,radii=NULL,std_meth="centers",
                        obs_names=NULL, var_names=NULL,
                        method=c("mrpca","cpca","vpca")){
 library(tidyverse)
  library(FactoMineR)
  source("int_centering.R") 
  source("radius_rotation.R")
  
  if(is.null(x)&&is.null(centers)){
    print("data must be provided either in intervals, or in centers and radii")
  }
  
  
  if(is.null(x)){
    C=centers
    R=radii
    
  }else{
    pp=ncol(x)
    x_min = x[,seq(from=1,to=pp,by=2)]
    x_max = x[,-seq(from=1,to=pp,by=2)]
    C=(x_min+x_max)/2
    R=(x_max-x_min)/2
  }
  
  if(is.null(var_names)){
    var_names=paste0("var_",1:ncol(C))
    obs_names=paste0("ob_",1:nrow(C))
    }

  
  out_int = int_centering(centers=data.matrix(C), radii=data.matrix(R), std_meth = std_meth)
  
  
  #############################
  ### SUPPLEMENTARY STUFF #####
  #############################
  #############################
  nsup<-0
  if (supp>0) {
    nsup<-supp 
    n<-dim(C)[1]-supp #se sono presenti punti supplementari si ridimensiona n
  }
  
  
  if (supp>0){
    Csupp<-C[(n+1):dim(C)[1],]
    Rsupp<-R[(n+1):dim(C)[1],]
    C<-C[1:n,]
    R<-R[1:n,]
    sup_out=int_centering(centers=Csupp, radii=Rsupp,std_meth = std_meth)
  
  }
  #############################
  #############################
  
  #################################
  #### DECOMPOSIZIONI       #######
  #### ( rivedere la notazione) ###
  #################################
  
  #############################
  #############################
  if(method=="mrpca"){
  cc=out_int$cc
  dd=out_int$dd
  
  cv2<-((abs(t(cc))%*%abs(dd))+(abs(t(dd))%*%abs(cc)))/nrow(C)
  ccc<-(t(cc)%*%cc)/n #matrice di varianza e covarianza dei centri finali
  inertia_centers<-sum(diag(ccc))
  ddd<-(t(dd)%*%dd)/n  #matrice di varianza e covarianza dei raggi finali 
  inertia_radii<-sum(diag(ddd))
  InTotale<-(ccc+ddd+cv2) #matrice di correlazione globale tra centri finali e raggi finali
  svd_ccc=svd(ccc)
  vv=svd_ccc$v
  v=svd_ccc$u
  d=svd_ccc$d
  svd_ddd=svd(ddd)
  uu<-svd_ddd$v
  b<-svd_ddd$d
  u<-svd_ddd$u
  autcen<-d
  autvetcen<-v
  autrag<-b
  autvetrag<-u
  nw1<-cc%*%v #proiezioni dei centri finali
  nw2<-dd%*%u #proiezione dei raggi finali
  CoCen<-nw1
  
  if(supp>0){
  ccsupp=sup_out$cc
  ddsupp=sup_out$dd
    nw1supp<-ccsupp%*%v 
    nw2supp<-ddsupp%*%u
  } #proiezione di centri e raggi finali per gli eventuali punti supplementari
  
  
  ragcord<-(t(dd)%*%dd%*%u%*%diag(b^-0.5))/n
  
  varcord<-t(cc)%*%cc%*%v%*%diag(d^-0.5)/n
  
    #################################
  #################################
  
  #################################
  ##### ROTAZIONI           #######
  #################################
  rot_out=radius_rotation(cc,dd)
  # rot_mat=rot_out$rotation_mat
  out = rot_out
  out$cen_coord=nw1
  out$no_rot_rad_coord=nw2
  
  ctr=sqrt(diag(out$last_rot_rad_coord[,1:2]%*%t(out$last_rot_rad_coord[,1:2])))/sqrt(diag(data.matrix(out$scaled_rad)%*%t(data.matrix(out$scaled_rad))))
  ctr=ctr/max(ctr)*1.01
  out$contributions = ctr
  out$center_correlations=varcord
  out$radii_correlations=ragcord
  out$variable_names=var_names
  out$observations_names=obs_names
  
  class(out)="interval_PCA_mrpca"
  }else if(method=="cpca"){
  out=list()
  
  pca_centers = PCA(cc) #starting analysis
  cc_svd=svd(cc)
  A=cc_svd$u %*% diag(cc_svd$d)
  B=cc_svd$v
  B_plus = B
  B_plus[B_plus<=0]=0
  B_minus = B
  B_minus[B_minus>=0]=0
  Z_plus=cc+dd
  Z_minus=cc-dd
  A_plus = Z_plus %*% B_minus + Z_minus %*% B_plus
  A_minus = Z_plus %*% B_plus + Z_minus %*% B_minus
  
  out$pca_centers = pca_centers
  out$A_plus=A_plus
  out$A_minus=A_minus
  
  class(out)="interval_PCA_cpca"
  
  }else if(method=="vpca"){
    
    out = list()
    print("coming soon")
    
    class(out)="interval_PCA_vpca"
  }else{
    print("please specify an implemented method")
  }
  
  
  return(out)
  
  
  
  #################################
  #################################
  ##### VISUALIZATION       #######
  #################################
  #################################
  
    
  }
