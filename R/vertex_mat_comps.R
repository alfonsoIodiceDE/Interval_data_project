rm(list=ls())
library("tidyverse")
source("interval_PCA.R")
source("plot.interval_PCA_cpca.R")
source("plot.interval_PCA_mrpca.R")
pepper = read_csv2(file="./Dataset/peperoni4vars.csv")



lab <- pepper[,1] #etichette delle osservazioni
pepper <- data.matrix(pepper)
pepper <- pepper[,-1]
interval_min <- pepper[,seq(1,7,2)] #matrice dei minimi
interval_max <- pepper[,seq(2,8,2)] #matrice dei massimi
C <- (interval_min+interval_max)/2 #matrice dei centri
R <- (interval_max-interval_min)/2 #matrice dei raggi



varnames <- c("H20","Protein","Lipid","Glucide") #nomi della variabili
interval_min=interval_min %>% as_tibble() 
names(interval_min) <-  varnames

interval_max=interval_max %>% as_tibble()
names(interval_max)<-varnames


source("vertex_mat_make.R")
interval_min = interval_min %>% mutate(what="min")
interval_max = interval_max %>% mutate(what="max")
out_vx=vertex_mat_make(interval_min = interval_min,interval_max = interval_max)


library("FactoMineR")
out=out_vx
pca_centers = PCA(out$vx_mat,scale.unit = FALSE, graph=FALSE)

out$pca_vx = pca_centers
out$inertias = pca_centers$eig
out$cos2s=pca_centers$ind$cos2
out$contr=pca_centers$ind$contr


out$obs_names=obs_names
out$var_names=var_names
