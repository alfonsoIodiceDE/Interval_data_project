rownames(Zc) <- lab
colnames(Zc) <- varnames

#PCA tradizionale sulla matrice Zc dei centri standardizzati
PCAc <- prcomp(Zc) 
#prcomp: effettua SVD 

library(factoextra)
fviz_eig(PCAc) #% di varianza spiegata da ogni dimensione 
eig.val <- get_eigenvalue(PCAc)
eig.val #autovalori e percentuali di varianza spiegata per ciascuna PC

#Valori rilevanti per le variabili
res.var <- get_pca_var(PCAc)
res.var$coord     #coordinate    
res.var$contrib   #contributi assoluti
res.var$cos2      #contributi relativi 

#Valori rilevanti per le osservazioni 
res.ind <- get_pca_ind(PCAc)
res.ind$coord     #coordinate    
res.ind$contrib   #contributi assoluti     
res.ind$cos2      #contributi relativi     

#Grafico delle unita' con contributi relativi = qualita' della rappresentazione sulle prime due dimensioni 
fviz_pca_ind(PCAc, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

#Grafico delle variabili con contributi assoluti alle prime due dimensioni 
fviz_pca_var(PCAc, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE )

#Biplot individui-variabili sulle prime due PC
fviz_pca_biplot(PCAc, repel = TRUE,
                col.var = "#2E9FDF", 
                col.ind = "#696969")

svdcc <- svd(Zc) #X = UDV'
u <- svdcc$u #autovettori di XX'
#In U sono contenuti gli score standardizzati delle PC, che rappresentano X nello spazio delle componenti principali
v <- svdcc$v #autovettori X'X
d <- svdcc$d #vettore dei valori singolari
all(round(u%*%diag(d)%*%t(v),6) == round(Zc,6)) #dimostrazione SVD

sigma <- diag(d)
B <- v #matrice dei loadings che esprimono la correlazione tra variabili originali e PC
A <- u%*%sigma #matrice degli scores: coordinate di ogni centro su ogni componente principale

#confronto di score e loadings appena calcolati con quelli forniti da prcomp
all(round(PCAc$x,5) == round(A,5))
all(round(PCAc$rotation,5) == round(B,5))

all(round(A%*%t(B),5) == round(Zc,5)) #approssima Zc (X=AB') 
all(round(Zc%*%B,5) == round(A,5)) #B e' ortonormale 

Bpos <- matrix(0,p,p) #Matrice dei loadings positivi
Bneg <- matrix(0,p,p) #Matrice dei loadings negativi

for (i in 1:nrow(B)){
  for(j in 1:ncol(B)){
    if (B[i,j]<0){
      Bneg[i,j] <- B[i,j]
    } 
  }
} 

for (i in 1:nrow(B)){
  for(j in 1:ncol(B)){
    if (B[i,j]>=0){
      Bpos[i,j] <- B[i,j]
    }
  }
}

stand_interval_min <- Zc - Zr #matrice dei minimi standardizzati
stand_interval_max <- Zc + Zr #matrice dei massimi standardizzati

#A- e A+ matrici dei limiti inferiore e superiore degli scores
amin <- stand_interval_max%*%Bneg + stand_interval_min%*%Bpos
amax <- stand_interval_max%*%Bpos + stand_interval_min%*%Bneg

library(ggplot2)
cpca_plot <- ggplot(data=NULL,
                    aes(x=res.ind$coord[,1], y=res.ind$coord[,2])) + 
  geom_point(colour="red") + scale_x_continuous(name="") +
  scale_y_continuous(name="") +
  geom_rect(data=NULL, mapping=aes(xmin=amin[,1], xmax=amax[,1], ymin=amin[,2], ymax=amax[,2]),
            fill='cyan4', color="black", alpha=0.4) +
  theme(legend.position="none",
        plot.title = element_text(face="bold", size = 14)) + 
  ggtitle("Centers Principal Component Analysis") + theme_bw() +
  geom_text(aes(label = lab))
cpca_plot