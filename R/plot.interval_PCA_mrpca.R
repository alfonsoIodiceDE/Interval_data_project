plot.interval_PCA_mrpca<-function(x,obs_names=NULL,var_names=NULL){

  library(ggrepel)
  library(ggforce)

  
  
  

    

  # graph 1: initial solution
  p <- ncol(x$first_rot_rad_coord)
  pp <- 2*p
  first_rot_rad_coord <- as_tibble(x$first_rot_rad_coord) 
  names(first_rot_rad_coord) <- paste0("fr_r_coord",1:p)
  
  cen_coord <- as_tibble(x$cen_coord*(-1)) 
  names(cen_coord) <-paste0("cen_coord",1:p)
  
  #first_plot
  
    initial_plot <- bind_cols(cen_coord,first_rot_rad_coord) %>% 
      pivot_longer(names_to = "elem", values_to="value", cols=1:all_of(pp)) %>% 
    separate(elem,sep="_c",into=c("what","dimension")) %>% 
    mutate(dimension=str_replace(dimension,pattern="oo",replacement = "coo"),) %>%
    group_by(what) %>% 
    pivot_wider(names_from = what,values_from=value,values_fn = list) %>% unnest(everything()) %>% 
    filter(dimension=="coord1"|dimension=="coord2") %>% 
    pivot_wider(names_from = dimension,values_from=2:3,values_fn = list) %>% unnest(everything()) %>% 
    mutate(obs_names=x$observations_names) %>% 
    ggplot(aes(x=cen_coord1,y=cen_coord2))+geom_point(col="blue",cex=1.5)+geom_text(aes(label=obs_names),cex=3.5)+
    geom_rect(aes(xmin=cen_coord1-fr_r_coord1,xmax=cen_coord1+fr_r_coord1,
                  ymin=cen_coord2-fr_r_coord2,ymax=cen_coord2+fr_r_coord2
    ),alpha=.25,fill="red",col="black")+theme_bw()+xlab("")+ylab("")+
    geom_segment(aes(x=cen_coord1,y=cen_coord2,xend=cen_coord1+fr_r_coord1,
                     yend=cen_coord2+fr_r_coord2),col="red")+ggtitle("Starting PCA map")
  
  
  
  last_rot_rad_coord <- as_tibble(x$last_rot_rad_coord) 
  names(last_rot_rad_coord) <- paste0("ls_r_coord",1:p)
  
  final_plot <- bind_cols(cen_coord,last_rot_rad_coord) %>% 
    pivot_longer(names_to = "elem",values_to="value",cols=1:ncol(.)) %>% 
    separate(elem,sep="_c",into=c("what","dimension")) %>% 
    mutate(dimension=str_replace(dimension,pattern="oo",replacement = "coo")) %>% 
    pivot_wider(names_from = what,values_from=value,values_fn = list) %>% unnest(everything()) %>% 
    filter(dimension=="coord1"|dimension=="coord2") %>% 
    pivot_wider(names_from = dimension,values_from=2:3,values_fn = list) %>% unnest(everything()) %>%
    mutate(contribution=x$contributions,
           obs_names=x$observations_names) %>% 
    ggplot(aes(x=cen_coord1,y=cen_coord2))+geom_point(col="blue",cex=1.5)+geom_text(aes(label=obs_names),cex=3.5)+
    geom_rect(aes(xmin=cen_coord1-ls_r_coord1,xmax=cen_coord1+ls_r_coord1,
                  ymin=cen_coord2-ls_r_coord2,ymax=cen_coord2+ls_r_coord2,
                  fill=contribution),alpha=.35,col="blue")+theme_bw()+xlab("")+ylab("")+
    geom_segment(aes(x=cen_coord1,y=cen_coord2,xend=cen_coord1+ls_r_coord1,
                     yend=cen_coord2+ls_r_coord2),col="blue") +ggtitle("Midpoint-radii Principal Component Analysis")
  
  
  
  ##### sovrapporre rotated vs non rotated ranges
  
  rad_coord_plot <- x$no_rot_rad_coord %>%as_tibble() %>% mutate(obs_names=x$observations_names) %>% 
    ggplot(aes(x=V1,y=V2))+geom_segment(aes(x=0,y=0,xend=V1,yend=V2),color="cyan4",lwd=0.8)+
    geom_text_repel(aes(label=obs_names),cex=3.5)+
    theme_bw()+ggtitle("Unrotated ranges coordinates")+xlab("")+ylab("")
  
  
  center_corr_plot <- x$center_correlations %>%as_tibble() %>% mutate(var_names=x$variable_names) %>% 
    ggplot(aes(x=V1,y=V2,label=var_names))+
    xlim(-1,1)+ylim(-1,1)+coord_equal()+ geom_text_repel()+
    geom_segment(aes(x=0,y=0,xend=V1,yend=V2),arrow = arrow(length=unit(.1,"inches"),type="closed"),color="cyan4")+
    geom_circle(aes(x0=0,y0=0,r=1),color="darkgrey",inherit.aes = FALSE)+theme_minimal()+
    geom_vline(aes(xintercept=0),color="grey")+geom_hline(aes(yintercept=0),color="grey")+
    ggtitle("centers correlations")+xlab("")+ylab("")
  
  
  radii_corr_plot <- x$radii_correlations %>%as_tibble() %>% mutate(var_names=x$variable_names) %>% 
    ggplot(aes(x=V1,y=V2,label=var_names))+ geom_text_repel()+
    xlim(-1,1)+ylim(-1,1)+coord_equal()+
    geom_segment(aes(x=0,y=0,xend=V1,yend=V2),arrow = arrow(length=unit(.1,"inches"),type="closed"),color="cyan4")+
    geom_circle(aes(x0=0,y0=0,r=1),color="darkgrey",inherit.aes = FALSE)+theme_minimal()+
    geom_vline(aes(xintercept=0),color="grey")+geom_hline(aes(yintercept=0),color="grey")+
    ggtitle("ranges correlations")+xlab("")+ylab("")
  
  out <- list()
  out$initial_plot <- initial_plot
  out$final_plot <- final_plot
  out$rad_coord_plot <- rad_coord_plot
  out$center_corr_plot <- center_corr_plot
  out$radii_corr_plot <- radii_corr_plot
  
  ### IJAS PAPER FIGURES not needed for the function
  # 
  # ggsave("MR_PCA_f_r_sol.pdf", final_plot, width=8, height=6)
  # ggsave("MR_PCA_rad_coords.pdf", rad_coord_plot, width=8, height=6)
  # ggsave("MR_PCA_cen_cor.pdf", center_corr_plot, width=6, height=6)
  # ggsave("MR_PCA_rad_cor.pdf", radii_corr_plot, width=6, height=6)
  
  return(out)
}
