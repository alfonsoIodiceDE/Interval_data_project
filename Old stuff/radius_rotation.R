radius_rotation<-function(scaled_C,scaled_R,first_only=F){
  library("tidyverse")  
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
  scale_C=cc
  scale_R=dd
  A<-scale_C/sqrt(n)
  B<-scale_R/sqrt(n)
  r<-dim(scale_C)[2] #numero di colonne della matrice dei centri 
  U<-matrix(0,r,r) #si crea una matrice (r x r) con elementi nulli 
  r1<-dim(scale_C)[1] #numero di righe della matrice dei centri
  
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
  u=svd(t(scale_R)%*% scale_R)$u
  
  Rnr<-scale_R%*%u%*%Ti #proiezioni dei raggi per matrice di rotazione iniziale Tc=K1%*%K2
  
  if(supp>0){
    Rnrsupp<-scale_Rsupp%*%u%*%Ti
  }
  
  
  
  W<-t(A)%*%B%*%(diag(diag(t(B)%*%B)^-0.5)) #esprime le covarianze centri-raggi normalizzate rispetto a B 
  #(matrice dei centri) che determina lo spazio in cui avviene la rotazione.
  AA<-t(A)%*%A #matrice di varianza e covarianza dei raggi standardizzati
  rho<-svd(AA)$d #valori singolari della matrice di varianza e covarianza dei raggi standardizzati
  rho<-max(rho) #il perno su cui avviare la rotazione ? l'autovalore massimo di A'A, in quanto gli 
  #autovalori sono funzione dell'inerzia totale.
  
  F_mat<-t(W)%*%Tc
  Tc[,diag(F_mat)<0]=Tc[,diag(F_mat)<0]*-1
  
  #se gli elementi della diagonale della matrice F_mat sono negativi, le colonne in Tc verranno 
  #moltiplicate per -1. In questo modo non si tiene conto del verso, perch? ci interessa solo 
  #l'angolo che non cambia se cambia il verso e conviene lavorare con tutti versi positivi. 
  
  
  
  fold<-0
  f<-sum(diag(t(W)%*%Tc)%*%diag(diag(t(Tc)%*%AA%*%Tc))^0.5)
  # eps<-.Machine$double.eps*100000
  eps<-.000000001
  counter<-0
  #se il miglioramento ? minore di eps mi fermo 
  #? stato deciso un valore fisso, ma sarebbe meglio considerare .Machine$double.eps per 10 o 100 o 1000
  #Se eps troppo grande l'algoritmo potrebbe fermarsi su punti piatti, quando dopo in realt? la funzione 
  #crescerebbe o decrescerebbe
  
  while((f>(fold+abs(f)*eps))&(counter<1000)){
    #delta<-f-fold
    counter<-counter+1
    
    pl_vec = diag(t(Tc)%*%AA%*%Tc)
    ql_vec = diag(t(W)%*%Tc)
    pl_ql_vec=(pl_vec)^(-3/2) * ql_vec 
    pl_ql_vec2=2*pl_vec^(-1/2)*(1/ql_vec)
    
    if(first_only==T){
      Tc_one=matrix(rep(Tc[,1],times=ncol(Tc)),ncol=ncol(Tc))
      myU_a=(AA%*%Tc- rho * Tc_one) %*% diag(pl_ql_vec)}
    else{ 
      myU_a=(AA%*%Tc- rho * Tc) %*% diag(pl_ql_vec) 
    }
    myU_b =  Tc %*% diag(diag(t(W) %*% W))%*% diag(pl_ql_vec2)
    myU_c = W %*% diag(pl_vec^(-1/2))
    U = myU_a - myU_b - myU_c
    U[,ql_vec==0]=0
    P<-svd(U)$u
    j<-svd(U)$d
    Q<-svd(U)$v
    T_mat<--P%*%Q #nel paper di Kier e Groenen T=-P%*%Q'
    #(pl^(-3/2))*(ql*(AA*Tc(:,l) - rho*Tc(:,1))) - (2*(pl^-0.5)*(1/ql)*(W(:,l)'*W(:,l)*Tc(:,l))) - (pl^-0.5)*W(:,l)
    F_mat<-t(W)%*%T_mat
    
    T_mat[ ,diag(F_mat)<0]=T_mat[ ,diag(F_mat)<0]*-1
    
    fold<-f
    f<-sum(diag(t(W)%*%T_mat)%*%diag(diag(t(T_mat)%*%AA%*%T_mat)^-0.5))
    OKT<-T_mat
    Tc<-T_mat
    
  }
  # sum(j) #somma degli autovalori di U
  # OKT
  # counter
  # 
  out=list()
  out$first_rotation_matrix = Ti
  out$last_rotation_matrix = OKT
  out$scaled_rad=scaled_R
  out$first_rot_rad_coord=Rnr
  out$last_rot_rad_coord = scale_R %*% u %*% OKT
  
  return(out)
}
