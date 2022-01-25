############################################################
####################### SIMULAZIONI #######################
############################################################
library(dplyr)
library(MASS)
set.seed(1234)
#m <- t(matrix(c( 3, -1.5,   0.5,  .5,
#                 -1.5,  3.0,   0,     0,
#                 -1.5, -1.5,  -0.5, -.5),4,3))

## matrice delle medie (centri)
m <- t(matrix(c(  0,     0,   0 ,  0,
                 -2,     2,   0,   0,
                  2,     -2 ,  0,   0),4,3))


s <- matrix(c(   1,  .4,   -.1,  .1, 
                .4,   1,     0,   0,
                -.1,  0,     1, -.7, 
                .1,   0,   -.7,   1),4,4) 

#matrice di rotazione var3 e var 4 neg correlati, le altre indipendenti (aumenta varianza var 3 e var4)
#due var inutili che si correlano tendono a mascherare la vera info effetto masking
#sia per centri che per raggi

s1 <- matrix(c(   1,  0,   .5,   0, 
                  0,  1,   0,   -.5,
                 .5,  0,     1,   0,
                 .0,  -.5,   0,   1),4,4) 
#solo per i raggi

I <- matrix(c(  1,  0,  0, 0, 
                0,  1,  0, 0,
                0,  0,  1, 0,
                0,  0,  0, 1),4,4)

one <- matrix(rep(1,6),6,1)

#alpha <-  .25   # quantità di Noise misurato in termini di percentuale della varianza .2 corrisponde al 20%
#tau   <- .8 # peso della variabilità dei centri se tau=0.8 allora 80% variabilità dai centri e 20% dai raggi

#alpha 0.2 0.5 0.8 
#tau 0.2 0.4 0.6
alpha<-c(0.2,0.5,0.8)
tau<-c(0.8,0.6,0.4)

simul<-list()
c<-0


for(i in 1:length(tau)){
  for(j in 1:length(alpha)){
    for(l in 1:100){
      for (k in 1:3) {
        if (k==1) {C   <- mvrnorm(6,m[k,],I,empirical = TRUE)
        col <- matrix(rep(1,6),6,1)} else {
          C <- rbind(C,mvrnorm(6,m[k,],I,empirical = TRUE))
          col <- rbind(col, matrix(rep(k,6),6,1)) } #tre medie diverse 
      }
      
      R <- matrix(runif(18*4,0,1),18,4)
      Rm <- apply(R,2,mean)
      Cm <- apply(C,2,mean)
      R <- R - matrix(1,18,1)%*%Rm
      C <- C - matrix(1,18,1)%*%Cm
      Y <- svd(R)
      Z <- svd(C)
      C <- Z$u%*%diag(Z$d) 
      R <- Y$u%*%diag(Y$d)
      D <- diag(runif(4,0,1))
      ones <- matrix(rep(1,18*4),18,4)
      
      
      #triangolarizzazione di chol
      F <- chol(s)
      B <- chol(s1)
      
      M  <- C%*%F
      S1 <- R%*%F
      S2 <- M%*%B + ones%*%D #per avere raggi correlati ai centri con rotazione rispetto a B
      
      Nm <- matrix(runif(18*4,0,2),18,4) #generare rumore
      Ns <- matrix(runif(18*4,0,.5),18,4)
      
      sNm <-(diag(diag(cov(Nm)))^.5)*tau[i]
      sNs <-(diag(diag(cov(Ns)))^.5)*(1-tau[i])
      
      M <- apply(M,2,scale)%*%sNm #attribuisce una porzione tau della varianza del rumore
      S1 <- apply(S1,2,scale)%*%sNs+matrix(rep(1,18),18,1)%*%Rm
      S2 <- apply(S2,2,scale)%*%sNs+matrix(rep(1,18),18,1)%*%Rm
      
      M <- M + alpha[j]*Nm #aggiunge rumore che sar? correlato a varianza del rumore
      
      S1 <- (S1 + alpha[j]*Ns)
      S2 <- (S2 + alpha[j]*Ns)
      
      for (r in 1:ncol(S1)){
        for (q in 1:nrow(S1)) {
          if (S1[q,r]<0) {S1[q,r]<-0}
          if (S2[q,r]<0) {S2[q,r]<-0}
        }} #elimina raggi negativi
      c<-c+1
      simul[[c]]<-list("alpha"=as.matrix(alpha[j]),"tau"=as.matrix(tau[i]),"M"=as.matrix(M),"S1"=as.matrix(S1),"S2"=as.matrix(S2))
    }
  }
}

######################################################################
############################ MRPCA FUNCTION ##########################
######################################################################

library(ggplot2)
library(ggforce)
library(ggrepel)

#### Inertial totale

mrpca<-function(C,R,lab,varnames){
  n <- dim(C)[1] #numero di osservazioni simboliche
  p <- dim(C)[2] #numero di variabili 
  Zc<-C
  Zr<-R
  cv2<-(abs(t(Zc)%*%Zr) + abs(t(Zr)%*%Zc))/n
  
  #cv2: componente della matrice dei correlazione globale relativa alle interconnessioni centri raggi
  covZc <- (t(Zc)%*%Zc)/n #matrice di varianza e covarianza dei centri
  InCentri <-sum(diag(covZc)) #inerzia dei centri
  covZr <- (t(Zr)%*%Zr)/n  #matrice di varianza e covarianza dei raggi
  InRaggi <- sum(diag(covZr)) #inerzia dei raggi
  InTotale <- (covZc + covZr + cv2) #matrice di correlazione globale 
  
  #InCentri/sum(diag(InTotale))
  #InRaggi/sum(diag(InTotale))
  
  #SVD di covZc e covZr per determinare le proiezioni di centri e raggi negli spazi generati dalle rispettive PC
  v <- svd(covZc)$u
  d <- svd(covZc)$d 
  b <- svd(covZr)$d 
  u <- svd(covZr)$u
  
  #Proiezioni di centri e raggi nei relativi spazi delle PC
  nw1 <- Zc%*%v 
  nw2 <- Zr%*%u 
  
  #Coordinate di centri e raggi
  ragcord <- (t(Zr)%*%Zr%*%u%*%diag(b^-0.5))/n
  varcord <- t(Zc)%*%Zc%*%v%*%diag(d^-0.5)/n
  
  # ----------------------------------------------------------------- #
  #  Determinazione della matrice di rotazione - INIZIO PROCEDURA --- #
  # ----------------------------------------------------------------- #
  #  Notazione Kiers e Groenen -------------------------------------- #
  # ----------------------------------------------------------------- #
  
  A <- Zc/sqrt(n)
  B <- Zr/sqrt(n)
  r <- dim(C)[2] 
  U <- matrix(0,r,r) 
  
  k <- svd(abs(t(A))%*%abs(B)+abs(t(B))%*%abs(A)) 
  #k: matrice di congruenza che mette in relazione blocchi di variabili diverse (centri e raggi) fornendo una sorta di covarianza centro-raggio
  K1 <- k$u
  K2 <- k$v
  Ti <- K1%*%K2  #prodotto tra le matrici dei vettori singolari di sinistra e di destra nella matrice k=(|A'||B|+|B'||A|)
  #K1 e K2 rappresentano delle direzioni, piu' sono simili tra loro piu' il loro prodotto dara' la matrice diagonale, mentre sarebbe una matrice piena se  fossero ortogonali (indipendenti). Si vuole massimizzare la traccia, in modo  da avere direzioni k1 e k2 simili per ogni variabile: si cercano connessioni forti tra il centro e il raggio di una variabile, non tra centri e raggi di variabili diverse. 
  Tc<-Ti
  #Tc e' l'inizializzazione ortonormale della matrice di rotazione T
  
  Rnr1<-Zr%*%u%*%Tc #proiezioni dei raggi per matrice di rotazione iniziale
  
  
  W <- t(A)%*%B%*%(diag(diag(t(B)%*%B)^-0.5)) #esprime le covarianze centri-raggi normalizzate rispetto alla matrice dei centri che determina lo spazio in cui  avviene la rotazione.
  AA <- t(A)%*%A #matrice var/cov dei raggi
  rho <- svd(AA)$d 
  rho <- max(rho) #il perno su cui avviare la rotazione e' l'autovalore massimo di A'A, in quanto gli autovalori sono funzione dell'inerzia totale.
  
  F <- t(W)%*%Tc
  for (i in 1:r){
    if (F[i,i]<0){
      Tc[,i] = Tc[,i]*-1
    } 
  }
  
  fold <- 0
  f <- sum(diag(t(W)%*%Tc)%*%diag(diag(t(Tc)%*%AA%*%Tc))^0.5)#criterio dell'algoritmo che si aggiorna ad ogni iterazione
  
  counter <- 0
  eps <- 0.000000001 #se il miglioramento e' minore di eps l'algoritmo convergera'
  
  
  while((f>(fold+abs(f)*eps))&(counter<1000)){
    counter<-counter+1
    
    pl_vec <- diag(t(Tc)%*%AA%*%Tc)
    ql_vec <- diag(t(W)%*%Tc)
    pl_ql_vec <-(pl_vec)^(-3/2) * ql_vec 
    pl_ql_vec2 <- 2*pl_vec^(-1/2)*(1/ql_vec)
    
    myU_a <- (AA%*%Tc- rho * Tc) %*% diag(pl_ql_vec) 
    myU_b <- Tc %*% diag(diag(t(W) %*% W))%*% diag(pl_ql_vec2)
    myU_c <- W %*% diag(pl_vec^(-1/2))
    U <- myU_a - myU_b - myU_c
    U[,ql_vec == 0] <- 0
    P <- svd(U)$u
    j <- svd(U)$d
    Q <- svd(U)$v
    T_mat <- -P%*%Q 
    F_mat <- t(W)%*%T_mat
    
    T_mat[ ,diag(F_mat)<0]<-T_mat[ ,diag(F_mat)<0]*-1
    
    fold <- f
    f <- sum(diag(t(W)%*%T_mat)%*%diag(diag(t(T_mat)%*%AA%*%T_mat)^-0.5))
    OKT <- T_mat
    Tc <- T_mat
  }
  
  OKT
  counter
  
  #Coordinate dei raggi ruotati in seguito alla rotazione procustiana
  Rnr2 <- Zr%*%OKT%*%u
  
  ctr <- sqrt(diag(Rnr2[,1:2]%*%t(Rnr2[,1:2])))/sqrt(diag(Zr%*%t(Zr))) 
  #ctr: contributi relativi sulle prime due dimensioni dati dal rapporto tra la lunghezza dei raggi ruotati e quella originale
  #i contributi relativi valutano la qualita' di rappresentazione delle unita'
  
  initial_plot <- ggplot(data=NULL, aes(x=nw1[,1], y=nw1[,2])) +
    geom_point(col="red", cex=2) +
    geom_text(aes(label=lab), cex=4, color="black") +
    geom_rect(aes(xmin=nw1[,1] - Rnr1[,1], xmax=nw1[,1] + Rnr1[,1],
                  ymin=nw1[,2] - Rnr1[,2], ymax=nw1[,2] + Rnr1[,2]),
              alpha=.3, fill="red", col="black") + theme_bw() + xlab("") +
    ylab("") + geom_segment(aes(x=nw1[,1], y=nw1[,2], xend=nw1[,1]+Rnr1[,1],
                                yend=nw1[,2]+Rnr1[,2]), col="darkblue")+
    theme(legend.position="none",
          plot.title = element_text(face="bold", size = 14))+
    ggtitle("Soluzione iniziale")
  initial_plot
  
  
  #ctr=ctr/max(ctr)*1.01
  contributions<-ctr
  
  final_plot <- ggplot(data=NULL, aes(x=nw1[,1], y=nw1[,2])) +
    geom_point(col="red", cex=2) +
    geom_text(aes(label=lab), cex=4, color="black")+
    geom_rect(aes(xmin=nw1[,1] - Rnr2[,1], xmax=nw1[,1] + Rnr2[,1], 
                  ymin=nw1[,2] - Rnr2[,2], ymax=nw1[,2] + Rnr2[,2],
                  fill=contributions), alpha=0.45, col="black") + 
    theme_bw() + xlab("") + ylab("") + 
    geom_segment(aes(x=nw1[,1], y=nw1[,2], xend=nw1[,1] + Rnr2[,1],
                     yend=nw1[,2] + Rnr2[,2]), col="darkblue")+
    theme(legend.justification="top",
          plot.title = element_text(face="bold", size = 14))+
    ggtitle("Soluzione finale")
  final_plot
  
  rad_coord_plot <- ggplot(data=NULL,aes(x=nw2[,1],y=nw2[,2])) +
    geom_segment(aes(x=0, y=0 ,xend=nw2[,1], yend=nw2[,2]), 
                 color="cyan4", lwd=1)+
    geom_text(aes(label=lab), cex=4) + xlab("") + ylab("") +
    theme_bw() + ggtitle("Coordinate dei raggi non ruotati")
  rad_coord_plot
  
  cor(Zr)
  cor(Rnr2)
  
  
  center_corr_plot <- ggplot(data=NULL, 
                             aes(x=varcord[,1], y=varcord[,2], label=varnames)) +
    xlim(-1,1) + ylim(-1,1) + coord_equal() + geom_text_repel(label=varnames) + xlab("") +
    ylab("") + geom_segment(aes(x=0, y=0, xend=varcord[,1], yend=varcord[,2]),
                            arrow = arrow(length=unit(.25,"cm"), type="open"), color="cyan4",
                            lwd=1) + geom_circle(aes(x0=0, y0=0, r=1), color="slategrey", 
                                                 inherit.aes = FALSE, lwd=1) + theme_minimal() +
    geom_vline(aes(xintercept=0), color="grey", lwd=0.8) + 
    geom_hline(aes(yintercept=0), color="grey")+
    ggtitle("Correlazioni dei centri")
  center_corr_plot
  
  radii_corr_plot <- ggplot(data=NULL, 
                            aes(x=ragcord[,1], y=ragcord[,2], label=varnames)) + 
    geom_text_repel(label=varnames) + xlim(-1,1) + ylim(-1,1) + coord_equal() + xlab("") +
    ylab("")+ geom_segment(aes(x=0, y=0, xend=ragcord[,1], yend=ragcord[,2]),
                           arrow = arrow(length=unit(.25,"cm"), type="open"), color="cyan4", lwd=1) +
    geom_circle(aes(x0=0, y0=0, r=1), color="slategrey", inherit.aes = FALSE, lwd=1) +
    theme_minimal()+ geom_vline(aes(xintercept=0), color="grey") + 
    geom_hline(aes(yintercept=0), color="grey")+
    ggtitle("Correlazioni dei raggi")
  radii_corr_plot
  
  globvar<-(t(Zc)%*%Zc + t(Zr)%*%Zr + abs(t(Zc)%*%Zr) + abs(t(Zr)%*%Zc))/n
  vars<-(t(nw1)%*%nw1+t(Rnr2)%*%Rnr2+abs(t(nw1)%*%Rnr2)+abs(t(Rnr2)%*%abs(nw1)))/n
  #matrice di varianza e covarianza tra le proiezioni dei centri nw1 e le proiezioni dei raggi ruotati 
  #secondo la rotazione procustiana 
  sum(diag(vars)) #inerzia spiegata vars
  resid<-sum(diag(globvar))-sum(diag(vars))
  In<-sum(diag(vars))*100/sum(diag(globvar))
  Inrag<-sum(diag(covZr))*100/sum(diag(globvar))
  
  
  out<-list("% Tot In"=In,"% Rad In"=Inrag,"iterations"=counter)
  return(out)
}



lab<-LETTERS[1:18]
varnames<-paste0("var",1:4)

C<-list()
R<-list()
mrpca.res <- matrix(NA,length(simul),5)

for (i in 1:length(simul)){
  C[[i]]<-simul[[i]][["M"]]
  R[[i]]<-simul[[i]][["S1"]] 
  mrpca.rest <- mrpca(C[[2]],R[[2]],lab,varnames)
  mrpca.res[i,1] <- mrpca.rest[1]
  mrpca.res[i,2] <- mrpca.rest[2]
  mrpca.res[i,3] <- mrpca.rest[3]
  mrpca.res[i,4] <- simul[[i]][1]
  mrpca.res[i,5] <- simul[[i]][1]
}  
  
  In_mean1<-mean(a1) #
  Rad_In_mean1<-mean(b1) #
  iterations_mean1<-mean(c1) #
  In_var1<-var(a1) #
  Rad_In_var1<-var(b1) #
  iterations_var1<-var(c1) #
}


for (i in 1:length(simul)){
  C[[i]]<-simul[[i]][["M"]]
  R[[i]]<-simul[[i]][["S2"]] 
  a2[i]<- unlist(mrpca(C[[i]],R[[i]],lab,varnames)[1]
  b2[i]<-unlist(mrpca(C[[i]],R[[i]],lab,varnames)[2])
  c2[i]<-unlist(mrpca(C[[i]],R[[i]],lab,varnames)[3])
  In_mean2<-mean(a2) #
  Rad_In_mean2<-mean(b2) #
  iterations_mean2<-mean(c2) #
  In_var2<-var(a2) #
  Rad_In_var2<-var(b2) #
  iterations_var2<-var(c2) #
}

for(i in 1:length(tau)){
  for(j in 1:length(alpha)){
    ma
    mb
    mc
    
