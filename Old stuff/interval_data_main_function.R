library(ggplot2)
library(R.matlab)

#C R,labels,varnames,supp

source("int_centering.R")
# if(is.null(supp)){
#   supp<-0
# }
C=data.matrix(C)
R=data.matrix(R)
n<-dim(C)[1] 
p<-dim(C)[2]
out=int_centering(centers=C, radii=R,method="centers")
Cc=out$Cc

Rc=out$Rc


#si procede allo stesso modo con eventuali punti supplementari
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
  sup_out=int_centering(centers=Csupp, radii=Rsupp,method="centers")
  Ccsupp=sup_out$Cc
  Cc1supp=sup_out$Cc1
  Rcsupp=sup_out$Rc
  Rc1supp=sup_out$Rc1
}

# total_covariance_matrix

vart<-(t(Cc)%*%Cc+t(Rc)%*%Rc+t(abs(Cc))%*%abs(Rc)+t(abs(Rc))%*%abs(Cc))/n
var1<-diag(vart)
cc=Cc/sqrt(var1)
dd=Rc/sqrt(var1)
oi=R/sqrt(var1)

if (supp>0) {
  ccsupp=Ccsupp/sqrt(var1)
  ddsupp=Rcsupp/sqrt(var1)
  oisupp=Rsupp/sqrt(var1)
}

# cv<-((abs(cc)%*%abs(t(dd)))+(abs(dd)%*%abs(t(cc)))) #? una matrice che fornisce le devianze e le 
#codevianze tra gli individui, che presenta la stessa traccia di cv2*n
cv2<-((abs(t(cc))%*%abs(dd))+(abs(t(dd))%*%abs(cc)))/nrow(C)

ccc<-(t(cc)%*%cc)/n #matrice di varianza e covarianza dei centri finali
inertia_centers<-sum(diag(ccc))
ddd<-(t(dd)%*%dd)/n  #matrice di varianza e covarianza dei raggi finali 
inertia_radii<-sum(diag(ddd))
InTotale<-(ccc+ddd+cv2) #matrice di correlazione globale tra centri finali e raggi finali

#Si effettua la SVD delle matrici var/cov dei centri finali e dei raggi finali e se ne 
#estraggono gli autovettori per le proiezioni nw1 e nw2.
svd_ccc=svd(ccc)
vv=svd_ccc$v
v=svd_ccc$u
d=svd_ccc$d

svd_ddd=svd(ddd)
uu<-svd_ddd$v
b<-svd_ddd$d
u<-svd_ddd$u
autcen<-d
autvetcen<-v
autrag<-b
autvetrag<-u

nw1<-cc%*%v #proiezioni dei centri finali
nw2<-dd%*%u #proiezione dei raggi finali
CoCen<-nw1

if(supp>0){
  nw1supp<-ccsupp%*%v 
  nw2supp<-ddsupp%*%u
} #proiezione di centri e raggi finali per gli eventuali punti supplementari

# ----------------------------------------------------------------- #
#  Determinazione della matrice di rotazione - INIZIO PROCEDURA --- #
# ----------------------------------------------------------------- #
#  Notazione Kiers e Groenen -------------------------------------- #
# ----------------------------------------------------------------- #
#Con la rotazione Procrustes una matrice A(n x r) viene rotata da una matrice ortonormale T (r x r)
#in modo da coincidere il pi? possibile con una matrice target B(n x r) in termini di minimi quadrati
#massimizzando il coefficiente di congruenza medio, che esprime il grado in cui mediamente le colonne sono
#proporziali. Questa ottimizzazione avviene attraverso la tecnica di majorization.
#Rotazione Procustiana: si effettua una rotazione dei raggi rispetto ai centri, massimizzando il 
#coefficiente di congruenza (Tucker 1951).


A<-cc/sqrt(n)
B<-dd/sqrt(n)
r<-dim(Cc)[2] #numero di colonne della matrice dei centri 
U<-matrix(0,r,r) #si crea una matrice (r x r) con elementi nulli 
r1<-dim(Cc)[1] #numero di righe della matrice dei centri

k<-svd(abs(t(A))%*%abs(B)+abs(t(B))%*%abs(A)) #matrice di congruenza che mette in relazione 
#blocchi di variabili diverse (centri e raggi) fornendo una sorta di covarianza centro-raggio
K1<-k$u
K2<-k$v
Ti<-K1%*%K2 #prodotto tra le matrici dei vettori singolari di sinistra e di destra nella matrice 
#(|A'||B|+|B'||A|).

#K1 e K2 rappresentano delle direzioni, pi? sono simili tra loro pi? il loro prodotto
#dar? la matrice diagonale, mentre sarebbe una matrice piena se fossero ortogonali (indipendenti).
#Si massimizzare la traccia, in modo da avere direzioni k1 e k2 simili per ogni variabile: si cercano
#connessioni forti tra il centro e il raggio di una variabile, non tra centri e raggi di variabili
#diverse. (non tra variabil diverse)
#e noi in effetti realizziamo una rotazione ortogonale
Tc<-Ti #Tc ? l'inizializzazione ortonormale della matrice di rotazione T


######### Determinazione di una soluzione iniziale##
#--------------------------------------------------#
#Determinazione delle coordinate dei raggi ruotati##
#--------------------------------------------------#
#################################################### 

Rnr<-dd%*%u%*%Ti #proiezioni dei raggi per matrice di rotazione iniziale Tc=K1%*%K2
if(supp>0){
  Rnrsupp<-ddsupp%*%u%*%Ti
}

x1min<-c()
y1min<-c()
x1max<-c()
y1max<-c()
for(i in 1:n){
  x1min[i]<-nw1[i,1]-Rnr[i,1]
  y1min[i]<-nw1[i,2]-Rnr[i,2]
  x1max[i]<-nw1[i,1]+Rnr[i,1]
  y1max[i]<-nw1[i,2]+Rnr[i,2]
}

X1<-c()
Y1<-c()
x1end<-c()
y1end<-c()
for (i in 1:n){
  X1[i]<-nw1[i,1]
  Y1[i]<-nw1[i,2]
  x1end[i]<-nw1[i,1]+Rnr[i,1]
  y1end[i]<-nw1[i,2]+Rnr[i,2]
}


initial_solution<-ggplot() + geom_point(mapping=aes(x=nw1[,1],y=nw1[,2]),colour="red")+
  scale_x_continuous(name="") + 
  scale_y_continuous(name="") +
  geom_rect(data=NULL, mapping=aes(xmin=x1min, xmax=x1max, ymin=y1min, ymax=y1max, fill="red"), color="black", alpha=0.5)+
  geom_segment(data=NULL,mapping=aes(x=X1,y=Y1,xend=x1end,yend=y1end),colour="blue")+
  theme(legend.position="none",plot.title = element_text(face="bold",size = 14))+
  ggtitle("Initial Solution")+
  geom_text(aes(x=nw1[,1],y=nw1[,2],label = lab),vjust="inward",hjust="inward")

initial_solution

W<-t(A)%*%B%*%(diag(diag(t(B)%*%B)^-0.5)) #esprime le covarianze centri-raggi normalizzate rispetto a B 
#(matrice dei centri) che determina lo spazio in cui avviene la rotazione.
AA<-t(A)%*%A #matrice di varianza e covarianza dei raggi standardizzati
rho<-svd(AA)$d #valori singolari della matrice di varianza e covarianza dei raggi standardizzati
rho<-max(rho) #il perno su cui avviare la rotazione ? l'autovalore massimo di A'A, in quanto gli 
#autovalori sono funzione dell'inerzia totale.

F<-t(W)%*%Tc
for (i in 1:r){
  if (F[i,i]<0){
    Tc[,i]=Tc[,i]*-1
  }#se gli elementi della diagonale della matrice F sono negativi, le colonne in Tc verranno 
  #moltiplicate per -1. In questo modo non si tiene conto del verso, perch? ci interessa solo 
  #l'angolo che non cambia se cambia il verso e conviene lavorare con tutti versi positivi. 
}

f<-0
fold<-0
f<-sum(diag(t(W)%*%Tc)%*%diag(diag(t(Tc)%*%AA%*%Tc))^0.5)

counter<-0
eps<-0.0000001 #se il miglioramento ? minore di eps mi fermo 
#? stato deciso un valore fisso, ma sarebbe meglio considerare .Machine$double.eps per 10 o 100 o 1000
#Se eps troppo grande l'algoritmo potrebbe fermarsi su punti piatti, quando dopo in realt? la funzione 
#crescerebbe o decrescerebbe
j<-0

while((f>(fold+abs(f)*eps))&(counter<1000)){
  #delta<-f-fold
  counter<-counter+1
  for (l in 1:r){
    pl<-t(Tc[,l])%*%AA%*%Tc[,l]
    ql<-t(W[,l])%*%Tc[,l]
    if(ql!=0){
      U[,l]<-as.numeric((pl^(-3/2)))*(as.numeric(ql)*(AA%*%matrix(Tc[,l]) - rho*matrix(Tc[,1]))) - 
        (2*as.numeric(pl^-0.5)*as.numeric(1/ql)*(as.numeric(t(matrix(W[,l]))%*%matrix(W[,l]))*matrix(Tc[,l]))) -
        as.numeric(pl^-0.5)*(W[,l])
    } else{
      U[,l]<-matrix(0,dim(U)[1],1)
    }
  }
  P<-svd(U)$u
  j<-svd(U)$d
  Q<-svd(U)$v
  T<--P%*%Q #nel paper di Kier e Groenen T=-P%*%Q'
  #(pl^(-3/2))*(ql*(AA*Tc(:,l) - rho*Tc(:,1))) - (2*(pl^-0.5)*(1/ql)*(W(:,l)'*W(:,l)*Tc(:,l))) - (pl^-0.5)*W(:,l)
  F<-t(W)%*%T
  for (i in 1:r){
    if(F[i,i]<0){
      T[,i]<-T[,i]*-1
    }
  }
  fold<-f
  f<-sum(diag(t(W)%*%T)%*%diag(diag(t(T)%*%AA%*%T)^-0.5))
  OKT<-T
  Tc<-T
  counter
}
sum(j) #somma degli autovalori di U
OKT
counter

#Determinazione delle coordinate dei raggi ruotati considerando la rotazione procustiana
Rnr<-dd%*%u%*%OKT
CoRag<-Rnr

if (supp>0){
  Rnrsupp<-ddsupp%*%u%*%OKT
}
##########################################SOLUZIONE FINALE##########################################

#colori<-abs(Rnr[,1])*abs(Rnr[,2])
#prodotto dei valori assoluti delle coordinate dei raggi ruotati 
#sulle prime due dimensioni
#colori<-colori/max(colori) #rapporto che indica 
colori<-sqrt(diag(Rnr[,1:2]%*%t(Rnr[,1:2])))/sqrt(diag(dd%*%t(dd)))
colori<-colori/(max(colori)*1.01)
#colori<-Rnr[,1:2]^2%*%t(Rnr[,1:2]^2)/dd[,1:2]^2%*%t(dd[,1:2]^2)
#misura la qualit? della rappresentazione: rapporto tra l'area dei rettangoli nel sottospazio e l'
#area dei rettangoli nello spazio originario

x2min<-c()
y2min<-c()
x2max<-c()
y2max<-c()
colore<-matrix(0,n,3)
for(i in 1:n){
  x2min[i]<-nw1[i,1]-Rnr[i,1]
  y2min[i]<-nw1[i,2]-Rnr[i,2]
  x2max[i]<-nw1[i,1]+Rnr[i,1]
  y2max[i]<-nw1[i,2]+Rnr[i,2]
  #colore[i,1]<-0.49
  #colore[i,2]<-1-colori[i]
  #colore[i,3]<-colori[i]
  colore[i,1]<-0.49
  colore[i,2]<-1-colori[i]
  colore[i,3]<-colori[i]
} #si disegnano i rettangoli i cui vertici sono rappresentati dagli estremi (min e max) sulle due 
#dimensioni

X2<-c()
Y2<-c()
x2end<-c()
y2end<-c()

for (i in 1:n){
  X2[i]<-nw1[i,1]  #coordinate proiezioni dei centri sulla prima dimensione
  Y2[i]<-nw1[i,2]  #coordinate proiezioni dei centri sulla seconda dimensione
  x2end[i]<-nw1[i,1]+Rnr[i,1] 
  y2end[i]<-nw1[i,2]+Rnr[i,2] 
}
#xend e  yend rappresentano punti dove il segmento che rappresenta il raggio dell'intervallo deve
#terminare sulla prima e sulla seconda dimensione, rispettivamente 
final_solution<-ggplot() + geom_point(mapping=aes(x=nw1[,1],y=nw1[,2]),colour="red")+
  scale_x_continuous(name="") + 
  scale_y_continuous(name="") +
  geom_rect(data=NULL, mapping=aes(xmin=x2min, xmax=x2max, ymin=y2min, ymax=y2max, fill=I(rgb(colore[,1], colore[,2],colore[,3]))), color="black", alpha=0.4)+
  geom_segment(data=NULL,mapping=aes(x=X2,y=Y2,xend=x2end,yend=y2end),colour="cyan4")+
  theme(legend.position="none",plot.title = element_text(face="bold",size = 14))+
  ggtitle("Final Solution")+
  geom_text(aes(x=nw1[,1],y=nw1[,2],label = lab),vjust="inward",hjust="inward")

xminsupp<-c()
yminsupp<-c()
xmaxsupp<-c()
ymaxsupp<-c()

if (supp>0){
  for (i in 1:nsupp){
    xminsupp[i]<-nw1supp[i,1]-Rnrsupp[i,1]
    yminsupp[i]<-nw1supp[i,2]-Rnrsupp[i,2]
    xmaxsupp[i]<-nw1supp[i,1]+Rnrsupp[i,1]
    ymaxsupp[i]<-nw1supp[i,2]+Rnrsupp[i,2]
  }
  Xsupp<-c()
  Ysupp<-c()
  xendsupp<-c()
  yendsupp<-c()
  
  for (i in 1:nsupp){
    Xsupp[i]<-nw1supp[i,1]
    Ysupp[i]<-nw1supp[i,2]
    xendsupp[i]<-nw1supp[i,1]+Rnrsupp[i,1]
    yendsupp[i]<-nw1supp[i,2]+Rnrsupp[i,2]
  }
  final_solution+geom_point(data=NULL,aes(nw1supp[,1],nw1supp[,2]),colour="blue")+
    geom_rect(data=NULL, mapping=aes(xmin=xminsupp, xmax=xmaxsupp, ymin=yminsupp, ymax=ymaxsupp, fill="yellow"), color="black", alpha=0.3)+
    geom_segment(data=NULL,mapping=aes(x=Xsupp,y=Ysupp,xend=xendsupp,yend=yendsupp),colour="green")+
    geom_text(aes(x=nw1[,1],y=nw1[,2],label = lab[n+1:dim(C)[1]]),vjust="inward",hjust="inward")
}
final_solution
#si ripete il ragionamento per gli eventuali punti supplementari
#########################################Coordinate dei raggi#########################################

X3<-c()
Y3<-c()
for (i in 1:n){
  X3[i]<-nw2[i,1]  #coordinate delle proiezioni dei raggi sulla prima dimensione
  Y3[i]<-nw2[i,2] #coordinate delle proiezioni dei raggi sulla seconda dimensione
}
radii_coord<-ggplot()+geom_segment(data=NULL,mapping=aes(x=0,y=0,xend=X3,yend=Y3),colour="cyan4",lwd=1.2)+
  scale_x_continuous(name="") + 
  scale_y_continuous(name="") +
  theme(legend.position="none",plot.title = element_text(face="bold",size = 14))+
  ggtitle("Radii coordinates")+
  geom_text(aes(x=nw2[,1],y=nw2[,2],label = lab))
radii_coord

if (supp>0){
  for (i in 1:nsupp){
    Xsupp[i]<-nw2supp[i,1]
    Ysupp[i]<-nw2supp[i,2]
  }
  return(radii_coord+geom_segment(data=NULL,mapping=aes(x=0,y=0,xend=Xsupp,yend=Ysupp),colour="green")+
           geom_text(aes(x=nw2supp[,1],y=nw2supp[,2],label = lab[n+1:dim(C)[1]]),vjust="inward",hjust="inward")) 
}
#si ripete il ragionamento per gli eventuali punti supplementari.
#######################################Calcolo dei residui#######################################
Rnr<-dd%*%u%*%OKT

if (supp>0){
  Rnrsupp<-ddsupp%*%u%*%OKT
}

vars<-(t(nw1)%*%nw1+t(Rnr)%*%Rnr+abs(t(nw1))%*%abs(Rnr)+abs(t(Rnr))%*%abs(nw1))/n
#matrice di varianza e covarianza tra le proiezioni dei centri nw1 e le proiezioni dei raggi ruotati 
#secondo la rotazione procustiana Rnr
sum(diag(vars)) #inerzia spiegata vars
resid<-dim(C)[2]-sum(diag(vars)) #differenza tra inerzia totale e quella della matrice vars

vars1<-(t(nw1[,1:2])%*%nw1[,1:2]+t(Rnr[,1:2])%*%Rnr[,1:2]+abs(t(nw1[,1:2]))%*%abs(Rnr[,1:2])+
          abs(t(Rnr[,1:2]))%*%abs(nw1[,1:2]))/n #si procede considerando solo le prime due dimensioni 
#per la possibilit? di rappresentazione grafica
In2CP<-sum(diag(vars1))*100/dim(C)[2] 
#proporzione di inerzia totale spiegata con due dimensioni

ragcord<--1*(t(dd)%*%dd%*%u%*%diag(b^-0.5))/n
ragcord[,1]^2+ragcord[,2]^2
varcord<--1*t(cc)%*%cc%*%v%*%diag(d^-0.5)/n

#######################################CORRELAZIONE DEI CENTRI#######################################

X4<-c()
Y4<-c()
for (i in 1:p){
  X4[i]<-varcord[i,1]
  Y4[i]<-varcord[i,2]
}
midpoints<-ggplot()+geom_segment(data=NULL,mapping=aes(x=0,y=0,xend=X4,yend=Y4),colour="cyan4",lwd=1.2)+
  scale_x_continuous(name="") + 
  scale_y_continuous(name="") +
  theme(legend.position="none",plot.title = element_text(face="bold",size = 14))+
  ggtitle("Correlation (Midpoints)")+
  geom_text(aes(x=varcord[,1],y=varcord[,2],label = varnames),size=6)

midpoints

#######################################CORRELAZIONE DEI RAGGI#######################################
X5<-c()
Y5<-c()
for (i in 1:p){
  X5[i]<-ragcord[i,1]
  Y5[i]<-ragcord[i,2]
}

radii<-ggplot()+geom_segment(data=NULL,mapping=aes(x=0,y=0,xend=X5,yend=Y5),colour="cyan4",lwd=1.2)+
  scale_x_continuous(name="")+
  scale_y_continuous(name="") +
  theme(legend.position="none",plot.title = element_text(face="bold",size = 14))+
  ggtitle("Correlation (Radii)")+
  geom_text(aes(x=ragcord[,1],y=ragcord[,2],label = varnames),size=6)

radii
res<-list(CoCen,CoRag,Ti,Tc,initial_solution,final_solution,radii_coord,midpoints,radii)
# return(res)



#########################################   PIPISTRELLI   #########################################
library("R.matlab")
dati<-readMat("C:\\Users\\schis\\Dropbox\\VivianaSchisa\\Rita\\batsall.mat")

C<-dati$dataCS #C ? la matrice dei centri standardizzati
R<-dati$dataRS #R ? la matrice dei raggi standardizzati
supp<-0 #nel dataset bats non sono presenti punti supplementari
labels<-dati$labels #sono i nomi delle 21 unit? statistiche
lab<-c()
for (i in 1:length(labels)){
  lab[i]<-labels[[i]]
  lab<-matrix(lab)
} #trasformazione della lista di liste "labels" in un vettore "lab" per semplicit?
idvar<-dati$idvar #sono i nomi delle quattro variabili: testa, coda, altezza e apertura alare
varnames<-c()
for (i in 1:length(idvar)){
  varnames[i]<-idvar[[i]]
  varnames<-matrix(varnames)
} #trasformazione della lista di liste "idvar" in un vettore "varnames" per semplicit?


##############################################   GRANO   ##############################################
mrs<-read.csv2("C:\\Users\\schis\\Desktop\\Interval Valued Data\\Grano - Centri e raggi standardizzati.csv",
               header=TRUE)
lab<-mrs[,1]
mrs<-data.matrix(mrs[,-1])
C<-mrs[,1:7]
R<-mrs[,8:14]
varnames<-c("humidity","weight","proteine","ash","gluten","iglut","yel")
supp<-0


real_input<-list(C,R,lab,varnames,supp)
Analisinew(real_input)

