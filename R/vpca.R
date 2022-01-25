#Si parte dalla matrice Z che alterna minimi e massimi standardizzati 

Zvert <- matrix(0,n*(2^p),p)	
for (i in 1:n){
  for (j in 1:p){
    for (k in 1:2^(p-j)){
      Zvert[((i-1)*2^p)+k,j]<-Z[i,2*j-1]
      Zvert[((i-1)*2^p)+2^(p-j)+k,j]<-Z[i,2*j]
    }
  }
}

for (l in 2:p){
  for (i in 1:n){
    for (j in l:p){
      for (k in 1:2^(p-j+l-1)){
        Zvert[((i-1)*2^p)+(2^(p-j+l-1))+k,j] <- Zvert[((i-1)*2^p)+k, j]
      }
    }
  }
}

svd <- svd(Zvert) #X=UDV'
u <- svd$u #autovettori XX'
v <- svd$v #autovettori X'X
d <- svd$d  #vettore dei valori singolari

sigma <- diag(d)
B <- v #matrice dei loadings
A <- u%*%sigma #matrice degli scores
all(round(u%*%diag(d)%*%t(v),5) == round(Zvert,5)) #dimostrazione SVD

colnames(Zvert)<-varnames
#PCA tradizionale sulla matrice Zc dei centri standardizzati
PCAv<-prcomp(Zvert) 
#prcomp: effettua i calcoli svd 

library(factoextra)

fviz_eig(PCAv) #% di varianza spiegata da ogni dimensione 
eig.valv <- get_eigenvalue(PCAv)
eig.valv #autovalori, % di varianza e % di varianza cumulata
#associati a ciascuna componente principale

#Valori rilevanti per le variabili
res.var <- get_pca_var(PCAv)
res.var$coord     #coordinate    
res.var$contrib   #contributi assoluti
res.var$cos2      #contributi relativi    

#Grafico delle variabili con contributi assoluti alle
#prime due dimensioni 
fviz_pca_var(PCAv, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
#Biplot individui-variabili sulle prime due PC
fviz_pca_biplot(PCAv,
                col.var = "#2E9FDF", 
                col.ind = "#696969") 
library(zoo)
q <- 2^p 
#Vertici minimi e massimi sulle prime due dimensionip per ogni unita'
xmin <- rollapply(A[,1],q,min,by=q)
xmax <- rollapply(A[,1],q,max,by=q)
ymin <- rollapply(A[,2],q,min,by=q)
ymax <- rollapply(A[,2],q,max,by=q)

cx <- (xmin+xmax)/2
cy <- (ymin+ymax)/2

library(ggplot2)

vpca_plot <- ggplot(data=NULL, aes(x=cx, y=cy)) + 
  geom_point(color="red", cex=3)+
  scale_x_continuous(name="") + 
  scale_y_continuous(name="") +
  geom_rect(data=NULL, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            fill="cyan4", color="black", alpha=0.4) +
  #geom_segment(aes(x=c1,y=c2,xend=rr1,yend=rr2),color="red") +
  theme(legend.position="none",
        plot.title = element_text(face="bold", size = 14)) +
  ggtitle("Vertex Principal Component Analysis") +
  geom_text(label = lab) + theme_bw()
vpca_plot