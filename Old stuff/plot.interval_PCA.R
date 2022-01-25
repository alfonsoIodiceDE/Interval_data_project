plot.interval_PCA<-function(x){
  library(tidyverse)
  library(ggrepel)
  library(ggforce)
  # x: interval PCA output object
  
  # out$first_rotation_matrix
  # out$last_rotation_matrix
  # out$first_rot_rad_coord
  # out$last_rot_rad_coord 
  # out$cen_coord
  # out$no_rot_rad_coord
  
  # graph 1: initial solution
  x=out
  p=ncol(x$first_rot_rad_coord)
  pp=2*p
  first_rot_rad_coord=as_tibble(x$first_rot_rad_coord) 
  names(first_rot_rad_coord)= paste0("fr_r_coord",1:p)
  
  cen_coord=as_tibble(x$cen_coord) 
  names(cen_coord)= paste0("cen_coord",1:p)
  
  #first_plot
  
  initial_plot=bind_cols(cen_coord,first_rot_rad_coord) %>% pivot_longer(names_to = "elem",values_to="value",cols=1:pp) %>% 
    separate(elem,sep="_c",into=c("what","dimension")) %>% 
    mutate(dimension=str_replace(dimension,pattern="oo",replacement = "coo"),) %>%
    group_by(what) %>% 
    pivot_wider(names_from = what,values_from=value,values_fn = list) %>% unnest(everything()) %>% 
    filter(dimension=="coord1"|dimension=="coord2") %>% 
    pivot_wider(names_from = dimension,values_from=2:3,values_fn = list) %>% unnest(everything()) %>% 
    mutate(obs_names=x$observations_names) %>% 
    ggplot(aes(x=cen_coord1,y=cen_coord2))+geom_point(col="blue",cex=.5)+geom_text(aes(label=obs_names),cex=2.5)+
    geom_rect(aes(xmin=cen_coord1-fr_r_coord1,xmax=cen_coord1+fr_r_coord1,
                  ymin=cen_coord2-fr_r_coord2,ymax=cen_coord2+fr_r_coord2
                  ),alpha=.25,fill="red",col="darkgrey")+theme_bw()+xlab("")+ylab("")+
    geom_segment(aes(x=cen_coord1,y=cen_coord2,xend=cen_coord1+fr_r_coord1,
                  yend=cen_coord2+fr_r_coord2),col="darkgrey")+ggtitle("interval PCA (starting) map")
  
  
  
  last_rot_rad_coord=as_tibble(x$last_rot_rad_coord) 
  names(last_rot_rad_coord)= paste0("ls_r_coord",1:p)
  
  
  
  final_plot=bind_cols(cen_coord,last_rot_rad_coord) %>% pivot_longer(names_to = "elem",values_to="value",cols=1:ncol(.)) %>% 
    separate(elem,sep="_c",into=c("what","dimension")) %>% 
    mutate(dimension=str_replace(dimension,pattern="oo",replacement = "coo")) %>% 
    pivot_wider(names_from = what,values_from=value,values_fn = list) %>% unnest(everything()) %>% 
    filter(dimension=="coord1"|dimension=="coord2") %>% 
    pivot_wider(names_from = dimension,values_from=2:3,values_fn = list) %>% unnest(everything()) %>%
    mutate(contribution=x$contributions,
           obs_names=x$observations_names) %>% 
    ggplot(aes(x=cen_coord1,y=cen_coord2))+geom_point(col="blue",cex=.5)+geom_text(aes(label=obs_names),cex=2.5)+
    geom_rect(aes(xmin=cen_coord1-ls_r_coord1,xmax=cen_coord1+ls_r_coord1,
                  ymin=cen_coord2-ls_r_coord2,ymax=cen_coord2+ls_r_coord2,
                  fill=contribution),alpha=.35,col="darkgrey")+theme_bw()+xlab("")+ylab("")+
    geom_segment(aes(x=cen_coord1,y=cen_coord2,xend=cen_coord1+ls_r_coord1,
                     yend=cen_coord2+ls_r_coord2),col="darkgrey") +ggtitle("interval PCA map (rotated)")
  
  
  
  
  rad_coord_plot=x$no_rot_rad_coord %>%as_tibble() %>% mutate(obs_names=x$observations_names) %>% 
    ggplot(aes(x=V1,y=V2))+geom_segment(aes(x=0,y=0,xend=V1,yend=V2),color="blue",alpha=.5)+
    geom_text(aes(label=obs_names),cex=2.5)+
    theme_bw()+ggtitle("(unrotated)radii coordinates map")
  
  
  center_corr_plot = x$center_correlations %>%as_tibble() %>% mutate(var_names=x$variable_names) %>% 
    ggplot(aes(x=V1,y=V2,label=var_names))+
    xlim(-1,1)+ylim(-1,1)+coord_equal()+ geom_text_repel()+
    geom_segment(aes(x=0,y=0,xend=V1,yend=V2),arrow = arrow(length=unit(.1,"inches"),type="closed"),alpha=.5)+
    geom_circle(aes(x0=0,y0=0,r=1),color="darkgrey",inherit.aes = FALSE)+theme_minimal()+
    geom_vline(aes(xintercept=0),color="grey")+geom_hline(aes(yintercept=0),color="grey")+ggtitle("centers correlations")
  
  
  radii_corr_plot = x$radii_correlations %>%as_tibble() %>% mutate(var_names=x$variable_names) %>% 
    ggplot(aes(x=V1,y=V2,label=var_names))+ geom_text_repel()+
    xlim(-1,1)+ylim(-1,1)+coord_equal()+
    geom_segment(aes(x=0,y=0,xend=V1,yend=V2),arrow = arrow(length=unit(.1,"inches"),type="closed"),alpha=.5)+
    geom_circle(aes(x0=0,y0=0,r=1),color="darkgrey",inherit.aes = FALSE)+theme_minimal()+
    geom_vline(aes(xintercept=0),color="grey")+geom_hline(aes(yintercept=0),color="grey")+ggtitle("radii correlations")
  
  # +
  #   xlab(paste0("d 1: ",round(inertia_table$perc_inertia[1],2),"%")) + 
  #   ylab(paste0("d 2: ",round(inertia_table$perc_inertia[2],2),"%"))
  # 
  
  
}
