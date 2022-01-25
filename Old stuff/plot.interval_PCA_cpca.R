plot.interval_PCA_cpca(x=NULL,ob_names=NULL,var_names=NULL){
  library(tidyverse)
  library(ggrepel)
  library(ggforce)
  library(janitor)
  # x: interval PCA output object
  
  if(is.null(var_names)){
    var_names <- paste0("var_",1:ncol(C))
  }
  
  if(is.null(obs_names)){
    obs_names <- paste0("ob_",1:nrow(C))  
  }
  
  
  
  x <- out
  center_coords = as_tibble(x$pca_centers$ind$coord) %>% clean_names()
  
  x$A_plus=x$A_plus %>% as_tibble() 
  names(x$A_plus)=paste0("dim_plus_",1:ncol(x$A_plus))
  
  x$A_minus=x$A_minus %>% as_tibble() 
  names(x$A_minus)=paste0("dim_minus_",1:ncol(x$A_minus))
  
  all_coords=cbind(center_coords,x$A_minus,x$A_plus)
 
  all_coords %>% ggplot(aes(x=dim_1,y=dim_2))+geom_point()+
  geom_rect(aes(xmin=dim_minus_1,xmax=dim_plus_1,ymin=dim_minus_2,ymax=dim_plus_2),alpha=.25)
   
}
