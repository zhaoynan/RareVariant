simulation_var<-function(n=n,cov=cov, E=E,G=G,EG=EG, res=res, model=model, mu=pis){
  nZ<-ncol(G)
  nX<-ncol(cov)+1
  if (model=='binomial') Vinv<-diag(mu*(1-mu))
  if (model=='gaussian') Vinv<-diag(rep(1,n))
  Xtu <- as.matrix(cbind(rep(1,n), cov, G, EG ))
  A <- t(Xtu) %*%  Vinv %*% Xtu
  D <- matrix(0,nrow=(nZ+nZ+nX),ncol=(nZ+nZ+nX))
  for (i in 1:n){
    Ui2<-as.vector(Xtu[i,] * res[i])
    D <- D+(Ui2)%*%t(Ui2)
  }
  A<-A/n
  A11<-A[1:(nX),1:(nX)]
  A12<-A[1:(nX),(nX+1):(nX+nZ+nZ)]
  A21<-A[(nX+1):(nX+nZ+nZ),1:(nX)]
  A22<-A[(nX+1):(nX+nZ+nZ),(nX+1):(nX+nZ+nZ)]
  D<-D/(n-1-nZ-nX)
  B11<-D[1:(nX),1:(nX)]
  B12<-D[1:(nX),(nX+1):(nX+nZ+nZ)]
  B21<-t(B12)
  B22<-D[(nX+1):(nX+nZ+nZ),(nX+1):(nX+nZ+nZ)]
  A11.edecomp <- eigen(A11,symmetric=T)  #general inverse
  A11.ev.inv <- ifelse(abs(A11.edecomp$values) >1e-8, 1/A11.edecomp$values, 0)
  A11inv <- t(A11.edecomp$vectors %*% (t(A11.edecomp$vectors) * A11.ev.inv))
  CovS <- (A22 - A21 %*% A11inv %*% A12 )
  CovS.edecomp <- eigen(CovS,symmetric=T)
  CovS.ev <- ifelse(abs(CovS.edecomp$values) >1e-8, CovS.edecomp$values, 0)
  V11 <- t(CovS.edecomp$vectors %*% (t(CovS.edecomp$vectors) * sqrt(CovS.ev)))
  return(V11)
}

calculation<-function(method=method, n.perm=n.perm, res=res, XdE=XdE, nZ=nZ, n=n, Us=Us, V11=V11, pow=pow, Ts_main=Ts_main,Ts_int=Ts_int,Ts=Ts, G=G, EG=EG){
  pPerm0_main = rep(NA,length(pow))
  pPerm0=rep(NA,length(pow))
  T0s_main = matrix(NA,n.perm,length(pow))
  P0s_main <-matrix(NA,n.perm,length(pow))
  pPerm0_int = rep(NA,length(pow))
  T0s_int = matrix(NA,n.perm,length(pow))
  P0s_int <-matrix(NA,n.perm,length(pow))
  T0s = matrix(NA,n.perm,length(pow))
  P0s <-matrix(NA,n.perm,length(pow))
  Y0 <- numeric(n)
  for (b in 1:n.perm){
    if (method=='ResidPerm') {res.p<-sample(res,n,replace=T)
    U0_main <- as.vector(t(G)%*% t(t(res.p)))
    U0_int <- as.vector(t(EG) %*% t(t(res.p)))
    } else if (method=='Bootstrap') {
      for(i in 1:n)  Y0[i] <- sample(c(1,0), 1, prob=c(pis[i], 1-pis[i]))
      tdat4=data.frame(trait=Y0,cov1)
      temp<-stats::glm(trait~.,data=tdat4,family=model)
      pis.p<-temp$fitted.value
      res.p<-tdat4$trait-pis.p
      U0_main <- as.vector(t(G)%*% t(t(res.p)))
      U0_int <-as.vector(t(EG) %*% t(t(res.p)))
    } else if (method=='Simulation'){
      nX <- ncol(XdE)+2
      Us <- rnorm(nZ+nZ)
      U0 <- (V11 %*% Us) * sqrt(n)
      U0_main <- U0[1:nZ]
      U0_int <- U0[(nZ+1):(nZ+nZ)]
    }
    
    for (j in 1:length(pow)){
      if (pow[j] < Inf){
        T0s_main[b,j] = round(sum((U0_main^pow[j])), digits = 8)
        T0s_int[b,j] = round(sum((U0_int^pow[j])), digits = 8)}
      if (pow[j] == Inf) {
        T0s_main[b,j] = round(max(abs(U0_main)), digits = 8)
        T0s_int[b,j] = round(max(abs(U0_int)), digits = 8)}
    } }
  
  T0s=T0s_main+T0s_int #unweighted version
  
  for (j in 1:length(pow)){
    pPerm0_main[j] = round((sum(abs(Ts_main[j])<=abs(T0s_main[1:(n.perm-1),j]))+1)/(n.perm), digits = 8) #why from B-1 and add 1-> to make p-value (0,1]
    P0s_main[,j] = (n.perm-rank(abs(T0s_main[,j]))+1)/(n.perm)
    
    pPerm0_int[j] = round((sum(abs(Ts_int[j])<=abs(T0s_int[1:(n.perm-1),j]))+1)/(n.perm), digits = 8)
    P0s_int[,j] = (n.perm-rank(abs(T0s_int[,j]))+1)/(n.perm)
    
    pPerm0[j] = round((sum(abs(Ts[j])<=abs(T0s[1:(n.perm-1),j]))+1)/(n.perm), digits = 8)
    P0s[,j] = (n.perm-rank(abs(T0s[,j]))+1)/(n.perm)
  }
  minp0_main<-apply(P0s_main,1,min)
  minp0_int<-apply(P0s_int,1,min)
  minp0<-apply(P0s,1,min)
  twomin<-cbind(minp0_main,minp0_int)
  min_min<-apply(twomin,1,min)
  fisher_min<-minp0_main*minp0_int
  
  Paspu_min<-(sum(min_min<=min(min(pPerm0_main),min(pPerm0_int)))+1)/(n.perm+1)
  Paspu_fisher<-(sum(log(fisher_min)<=log(min(pPerm0_main))+log(min(pPerm0_int)))+1)/(n.perm+1) #update log here for fisher apr12
  Paspu_unweighted<-(sum(minp0<=min(pPerm0))+1)/(n.perm+1)
  
  Ts <- c(Ts,min(pPerm0),min(min(pPerm0_main),min(pPerm0_int)),min(pPerm0_main)*min(pPerm0_int) )
  pvs <- c(pPerm0,Paspu_unweighted, Paspu_min, Paspu_fisher)
  names(Ts) <- c(paste("aGEj", pow, sep = ""), "aGEj","aGEj_minp","aGEj_fisher")
  names(pvs) = names(Ts)
  re <- list(min_min=min_min,fisher_min=fisher_min,pPerm0_main=pPerm0_main,pPerm0_int=pPerm0_int,  pvs=pvs)
  return(re)
}

