interval_PCA<- function(x=NULL,centers=NULL,radii=NULL,std_meth="centers",
                        obs_names=NULL, var_names=NULL,
                        method=c("mrpca","cpca","vpca"),supp=0){
  
  source("R/int_standardization.R")
  source("R/radius_rotation.R")
  source("R/vertex_mat_make.R")
  source("R/plot.interval_PCA_cpca.R")
  source("R/plot.interval_PCA_vpca.R")
  source("R/plot.interval_PCA_mrpca.R")
  
  
  # n <- dim(centers)[1] 
  # p <- dim(centers)[2]
  
  
  if(is.null(x)&&is.null(centers)){
    print("data must be provided either in intervals, or in centers and radii")
  }
  
  if(method!="vpca"){ #not vpca: need C and R
    if(is.null(x) ){ # x null -> C and R are provided already
      C=centers
      R=radii
      if(is.null(var_names)){
        var_names=paste0("var_",1:ncol(C))
        obs_names=paste0("ob_",1:nrow(C))
      }
    }else{ # x not null, not vpca: take minima and maxima from x, and compute C and R 
      
      pp=ncol(x)
      x_min = x %>% dplyr::select(seq(from=1,to=pp,by=2)) %>% as.matrix()
      x_max = x %>% dplyr::select(seq(from=2,to=pp,by=2)) %>% as.matrix()
      C=(x_min+x_max)/2
      R=(x_max-x_min)/2
      if(is.null(var_names)){
        var_names=paste0("var_",1:ncol(C))
        obs_names=paste0("ob_",1:nrow(C))
      }
    }
  }else{ # vpca case
    if(is.null(x) ){ #you need x, if it's null, compute it from C and R
      if(is.null(var_names)){
        var_names=paste0("var_",1:ncol(C))
        obs_names=paste0("ob_",1:nrow(C))
      }
      interval_min = as_tibble(C-R) %>% mutate(what="min")
      interval_max = as_tibble(C+R) %>% mutate(what="max")
      # x=rbind(interval_min,interval_max) %>%
      #   pivot_longer(cols=1:4,names_to = "var_names",values_to="value") %>% 
      #   arrange(var_names) %>% unite("var_names",what:var_names,sep = "_",remove = TRUE) %>% 
      #   pivot_wider(names_from=var_names,values_from = value,values_fn = list) %>% unnest(cols=everything())
      
    }else{ 
      if(is.null(var_names)){
        var_names=paste0("var_",1:(ncol(x)/2))
        obs_names=paste0("ob_",1:(nrow(x)))
      }# you have x
      pp=ncol(x)
      x_min = x %>% dplyr::select(seq(from=1,to=pp,by=2)) %>%as_tibble() 
      x_max = x %>% dplyr::select(seq(from=2,to=pp,by=2)) %>%as_tibble() 
      C = (x_min + x_max)/2
      R = (x_max - x_min)/2
    }
  }    
  
  out_int = int_standardization(x=NULL,centers=data.matrix(C), radii=data.matrix(R), std_meth = std_meth)

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
    sup_out=int_standardization(centers=Csupp, radii=Rsupp,std_meth = std_meth)
    
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
    
    n<-nrow(cc)
    p<-ncol(cc)
    
    cv2<-((abs(t(cc))%*%abs(dd))+(abs(t(dd))%*%abs(cc)))/n
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
  
    nw1<-cc%*%v #proiezioni dei centri finali
    nw2<-dd%*%u #proiezione dei raggi finali

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
    
    rot_out=radius_rotation(scaled_C=cc,scaled_R=dd)
    vars<-(t(nw1)%*%nw1+t(rot_out$last_rot_rad_coord)%*%rot_out$last_rot_rad_coord+abs(t(nw1)%*%rot_out$last_rot_rad_coord)+abs(t(rot_out$last_rot_rad_coord)%*%abs(nw1)))/n
    In<-sum(diag(vars))*100/sum(diag(out_int$glob_var))
    
    out = rot_out
    out$cen_coord=nw1
    out$no_rot_rad_coord=nw2
    out$mrpca_inertia=In
    
    ctr=sqrt(diag(out$last_rot_rad_coord[,1:2]%*%t(out$last_rot_rad_coord[,1:2])))/sqrt(diag(data.matrix(out$scaled_rad)%*%t(data.matrix(out$scaled_rad))))
    ctr=ctr/max(ctr)*1.01
    
    out$contributions = ctr
    out$center_correlations=varcord
    out$radii_correlations=ragcord
    out$variable_names=var_names
    out$observations_names=obs_names
    out$inertia_centers=inertia_centers
    out$inertia_radii=inertia_radii
    out$scaled_C=out_int$cc
    out$scaled_R=out_int$dd
    out$corr=out_int$Corr
    
    class(out)="interval_PCA_mrpca"
  }
  else if(method=="cpca"){
    out=list()
    cc=out_int$cc
    dd=out_int$dd
    
    pca_centers = PCA(cc,scale.unit = FALSE, graph=FALSE) #starting analysis
    
    
    cc_svd=svd(cc)
    A=cc_svd$u %*% diag(cc_svd$d)
    
    signs_check=colSums(abs((A[,1:2]>0)-(pca_centers$ind$coord[,1:2]>0)))
    A[,signs_check!=0]=A[,signs_check!=0]*(-1)
    
    
    B=cc_svd$v
    B[,signs_check!=0]=B[,signs_check!=0]*(-1)
    
    
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
    out$inertias = pca_centers$eig
    out$cos2s=pca_centers$ind$cos2
    out$contr=pca_centers$ind$contr
    out$obs_names=obs_names
    out$var_names=var_names
    out$scaled_C=out_int$cc
    out$scaled_R=out_int$dd
    out$corr=out_int$Corr
    class(out)="interval_PCA_cpca"
    
  }else if(method=="vpca"){
    
    
    pp=ncol(x)
    interval_min = as_tibble(out_int$cc-out_int$dd) %>% mutate(what="min")
    interval_max = as_tibble(out_int$cc+out_int$dd) %>% mutate(what="max")
    names(interval_max)=names(interval_min)
    vx_out = vertex_mat_make(interval_min = interval_min,interval_max = interval_max)
    vx_mat = vx_out$vx_mat
    vx_id=vx_out$ids
    
    
    out=list()
    pca_centers = PCA(vx_mat,scale.unit = FALSE, graph=FALSE)
    
    out$pca_vx = pca_centers
    out$inertias = pca_centers$eig
    out$cos2s=pca_centers$ind$cos2
    out$contr=pca_centers$ind$contr
    out$ids=vx_id
    out$obs_names=obs_names
    out$var_names=var_names
    out$scaled_C=out_int$cc
    out$scaled_R=out_int$dd
    out$corr=out_int$Corr
    
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
