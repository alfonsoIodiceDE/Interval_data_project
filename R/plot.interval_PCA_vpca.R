plot.interval_PCA_vpca<-function(x=NULL){
  
    
  mycoords = x$pca_vx$ind$coord %>%as_tibble %>% dplyr::select(Dim.1,Dim.2) %>%
    mutate(id=as_factor(x$ids)) %>% group_by(id) %>% 
    summarise(x_min = min(Dim.1),
              x_max = max(Dim.1),
              y_min = min(Dim.2),
              y_max = max(Dim.2)) %>% mutate(obs_names=x$obs_names)
  
  # vpca_plot=mycoords %>% 
  #   ggplot() + geom_rect(aes(xmin=x_min,xmax=x_max,ymin=y_min,ymax=y_max), 
  #                          fill="cyan4", color="black", alpha=0.4)
  
  
  vpca_plot = mycoords %>% ggplot(aes(x=(x_min+x_max)/2,y=(y_min+y_max)/2)) + geom_text(aes(label=obs_names),cex=3.5)+
    geom_rect(aes(xmin=x_min,xmax=x_max,ymin=y_min,ymax=y_max),fill="cyan4",alpha=.35,col="blue")+
    theme_bw() + ggtitle("Vertices Principal Component Analysis")+theme(legend.position="none") + 
    xlab("Dim_1") + ylab("Dim_2")
  
  out=list()
  out$coords=mycoords
  out$vpca_plot = vpca_plot

  return(out)
}
