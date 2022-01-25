C <- data.matrix(C)
R <- data.matrix(R)
n <- dim(C)[1] 
p <- dim(C)[2]

out_stand <- int_standardization(centers=C, radii=R, std_meth="symb")
cc <- out_stand$cc
dd <- out_stand$dd
Corr <- out_stand$Corr
glob_var <- out_stand$glob_var

#si procede allo stesso modo con eventuali punti supplementari
nsup <- 0
if (supp>0) {
  nsup <- supp 
  n <- dim(C)[1]-supp #se sono presenti punti supplementari si ridimensiona n
}


if (supp>0){
  Csupp <- C[(n+1):dim(C)[1],]
  Rsupp <- R[(n+1):dim(C)[1],]
  C <- C[1:n,]
  R <- R[1:n,]
  sup_out <- int_centering(x=NULL,centers=Csupp, radii=Rsupp,method="centers")
  Ccsupp <- sup_out$Cc
  Cc1supp <- sup_out$Cc1
  Rcsupp <- sup_out$Rc
  Rc1supp <- sup_out$Rc1
}

out_rot <- radius_rotation(scaled_C=cc,scaled_R=dd,first_only = F)


out_int_pca <- interval_PCA(centers=C,radii=R,std_meth="symb",method='mrpca')
varnames <- c("H20","Protein","Lipid","Glucide") 
mrpca_plot <- plot(out_int_pca,varnames)

out_cpca <- interval_PCA(centers=C,radii=R,std_meth="symb",method='cpca')
cpca_plot<-plot(out_cpca)
