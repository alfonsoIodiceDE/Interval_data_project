interval_data_maker <- function(x){
  library(MASS)
  library(tidyverse)
  alpha = x[1] %>% as.numeric 
  tau = x[2] %>% as.numeric 
  l = x[3] %>% as.numeric 

  # print("alpha,tau, rep")
  # print(c(alpha,tau, l))
  ############################################################
  ############################################################
  ## INITIALIZATION   ########################################
  ############################################################
  ############################################################
  ## matrice delle medie (centri)
  
  m <- t(matrix(c(  0,     -1,   0 ,  0,
                    -2,     2,   0,   0,
                    2,     2 ,  0,   0),4,3))
  
  
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
  
  for (k in 1:3) {
    
    if (k==1) {
      C   <- mvrnorm(6,m[k,],I,empirical = TRUE)
      col <- matrix(rep(1,6),6,1)
    }else{
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
Fm <- chol(s)
B <- chol(s1)

M  <- C%*%Fm
S1 <- R%*%Fm
S2 <- M%*%B + ones%*%D #per avere raggi correlati ai centri con rotazione rispetto a B

Nm <- matrix(runif(18*4,0,2),18,4) #generare rumore
Ns <- matrix(runif(18*4,0,.5),18,4)

sNm <-(diag(diag(cov(Nm)))^.5)*tau
sNs <-(diag(diag(cov(Ns)))^.5)*(1-tau)

M <- apply(M,2,scale)%*%sNm #attribuisce una porzione tau della varianza del rumore
S1 <- apply(S1,2,scale)%*%sNs+matrix(rep(1,18),18,1)%*%Rm
S2 <- apply(S2,2,scale)%*%sNs+matrix(rep(1,18),18,1)%*%Rm

M <- M + alpha*Nm #aggiunge rumore che sar? correlato a varianza del rumore

S1 <- (S1 + alpha*Ns)
S2 <- (S2 + alpha*Ns)

for (r in 1:ncol(S1)){
  for (q in 1:nrow(S1)) {
    if (S1[q,r]<0) {S1[q,r]<-0}
    if (S2[q,r]<0) {S2[q,r]<-0}
  }} #elimina raggi negativi

out=list()
out$M=M
out$S1=S1
out$S2=S2
return(out)
}
