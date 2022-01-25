plot.interval_PCA_cpca<-function(x=NULL){
  library(tidyverse)
  library(ggrepel)
  library(ggforce)
  library(janitor)
  # x: interval PCA output object
  
  
  
  
  center_coords = as_tibble(x$pca_centers$ind$coord) %>% clean_names()
  
  x$A_plus=x$A_plus %>% as_tibble() 
  names(x$A_plus)=paste0("dim_plus_",1:ncol(x$A_plus))
  
  x$A_minus=x$A_minus %>% as_tibble() 
  names(x$A_minus)=paste0("dim_minus_",1:ncol(x$A_minus))
  
  tot_contr = rowSums(x$contr[,1:2])
  tot_cos2 = rowSums(x$cos2s[,1:2])
  all_coords=cbind(center_coords,x$A_minus,x$A_plus,obs_names=x$obs_names) %>% #%>% pull(Id)
    mutate(
      contributions = tot_contr,
      cosine_sq = tot_cos2
    )
  
  centers_plot= all_coords %>% ggplot(aes(x=dim_1,y=dim_2)) + geom_text(aes(label=obs_names),cex=3.5,hjust=0,vjust=0) + 
    geom_point(aes(color=tot_cos2),cex=2)+ #size=contributions
    theme_bw()+geom_vline(xintercept=0,color="grey")+geom_hline(yintercept=0,color="grey") + 
    ggtitle("Centers Scatter Plot") + xlab("Dim_1") + ylab("Dim_2")
  
  attcoord = tibble(d_1=x$pca_centers$var$coord[,1],
                    d_2=x$pca_centers$var$coord[,2],
                    attribute=x$var_names,
                    ctr=rowSums(x$pca_centers$var$contrib[,1:2]),
                    cos2=rowSums(x$pca_centers$var$cos2[,1:2])
                    )
  
  corr_circle_plot = attcoord %>% ggplot(aes(x=d_1,y=d_2,label=attribute))+
    geom_text_repel()+xlim(-1,1)+ylim(-1,1)+coord_equal()+
    geom_segment(aes(x=0,y=0,xend=d_1,yend=d_2,color=ctr), #,alpha=cos2
                 arrow = arrow(length=unit(.1,"inches"),type="closed"))+
    geom_circle(aes(x0=0,y0=0,r=1),color="darkgrey",inherit.aes = FALSE)+theme_minimal()+
    geom_vline(aes(xintercept=0),color="grey")+geom_hline(aes(yintercept=0),color="grey") + 
    ggtitle("Loading plot") + xlab("Dim_1") + ylab("Dim_2")
  
  # +
  #   xlab(paste0("d 1: ",round(inertia_table$perc_inertia[1],2),"%")) + 
  #   ylab(paste0("d 2: ",round(inertia_table$perc_inertia[2],2),"%"))
  # 
  
  rect_plot=all_coords %>% ggplot(aes(x=dim_1,y=dim_2)) + geom_point(col="red",cex=1.5) +
    geom_text(aes(label=obs_names),cex=3.5)+
  geom_rect(aes(xmin=dim_minus_1,xmax=dim_plus_1,ymin=dim_minus_2,ymax=dim_plus_2),fill="cyan4",alpha=.35,col="blue")+
    theme_bw() + ggtitle("Centers Principal Component Analysis")+theme(legend.position="none") + 
    xlab("Dim_1") + ylab("Dim_2")
  
  out=list()
  
  out$centers_plot = centers_plot
  out$corr_circle_plot=corr_circle_plot
  out$rect_plot = rect_plot
  
  return(out)
}
