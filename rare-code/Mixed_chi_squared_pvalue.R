Get_Liu_PVal_Mod <- function(Q, W){
  Q <- as.matrix(Q)
  W <- as.matrix(W)
  
  A1 <- W
  A2 <- A1 %*% A1
  c1<-rep(0,4)
  c1[1] <- sum(diag(A1))
  c1[2] <- sum(diag(A2))
  c1[3] <- sum(A1*t(A2))
  c1[4] <- sum(A2*t(A2))
  param <- Get_Liu_Params_Mod(c1)
  
  Q.Norm <- (Q - param$muQ)/param$sigmaQ
  Q.Norm1 <- Q.Norm * param$sigmaX + param$muX
  p.value <- pchisq(Q.Norm1,df = param$l,ncp=param$d, lower.tail=FALSE)
  
  return(list(p.value = p.value[1], param=param,Q.Norm=Q.Norm1))
}

Get_Liu_Params_Mod <- function(c1){
  ## Helper function for getting the parameters for the null approximation
  
  muQ <- c1[1]
  sigmaQ <- sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2
  
  beta1<-sqrt(8)*s1
  beta2 <- 12*s2
  type1 <- 0
  
  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1 <- 1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <- l+d
  sigmaX <- sqrt(2)*a
  
  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}

Get_Matrix_Square.1 <- function(A){
  
  out <- eigen(A,symmetric=TRUE)
  ID1 <- which(out$values > 0)
  if(length(ID1)== 0){
    stop("Error to obtain matrix square!")
  }
  out1 <- t(out$vectors[,ID1]) * sqrt(out$values[ID1])
  return(out1)
}

Get_Davies_PVal <- function(Q, W){
  ##added by Zhangchen, for sparse matrix, 12.17.2018
  Q <- as.matrix(Q)
  W <- as.matrix(W)
  
  re <- Get_PValue(W, Q)
  param <- list()
  param$liu_pval<-re$p.val.liu[1]
  param$Is_Converged <- re$is_converge[1]
  
  re <- list(p.value = re$p.value[1], param=param, pval.zero.msg=re$pval.zero.msg )  
  return(re)
}

Get_PValue <- function(W,Q){
  
  lambda <- Get_Lambda(W)
  re <- Get_PValue.Lambda(lambda, Q)
  return(re)
}

Get_Lambda <- function(W){
  
  out.s <- eigen(W,symmetric=TRUE, only.values = TRUE)
  lambda1 <- out.s$values
  IDX1 <- which(lambda1 >= 0)
  IDX2 <- which(lambda1 > mean(lambda1[IDX1])/100000)
  
  if(length(IDX2) == 0){
    stop("No Eigenvalue is bigger than 0!!")
  }
  lambda <- lambda1[IDX2]
  return(lambda)
}

Get_PValue.Lambda <- function(lambda, Q){
  
  n1 <- length(Q)
  
  p.val <- rep(0,n1)
  Q.Norm <- rep(0,n1)
  p.val.liu <- rep(0,n1)
  is_converge <- rep(0,n1)
  p.val.liu <- Get_Liu_PVal.MOD.Lambda(Q, lambda)$p.value
  Q.Norm.liu <- Get_Liu_PVal.MOD.Lambda(Q, lambda)$Q.Norm
  for(i in 1:n1){
    out <- SKAT_davies(Q[i], lambda, acc=10^(-6))
    
    p.val[i] <- out$Qq
    is_converge[i]<-1
    
    # check convergence
    if(length(lambda) == 1){
      p.val[i] <- p.val.liu[i]
      Q.Norm[i] <- Q.Norm.liu[i]
    } else if(out$ifault != 0){
      is_converge[i]<-0
    }
    
    # check p-value
    if(p.val[i] > 1 || p.val[i] <= 0 ){
      is_converge[i]<-0
      p.val[i] <- p.val.liu[i]
      Q.Norm[i] <- Q.Norm.liu[i]
    }
  }
  
  p.val.msg = NULL
  p.val.log=NULL
  #cat(p.val[1])
  if(p.val[1] == 0){
    
    param <- Get_Liu_Params_Mod_Lambda(lambda)
    p.val.msg <- Get_Liu_PVal.MOD.Lambda.Zero(Q[1], param$muQ, param$muX, param$sigmaQ, param$sigmaX, param$l, param$d)
    p.val.log <- Get_Liu_PVal.MOD.Lambda(Q[1], lambda, log.p=TRUE)$p.value[1]
    
  }
  
  return(list(p.value=p.val, p.val.liu=p.val.liu, is_converge=is_converge, p.val.log=p.val.log, pval.zero.msg=p.val.msg))
  
}

Get_Liu_PVal.MOD.Lambda <- function(Q.all, lambda, log.p=FALSE){
  
  param <- Get_Liu_Params_Mod_Lambda(lambda)
  
  Q.Norm <- (Q.all - param$muQ)/param$sigmaQ
  Q.Norm1 <- Q.Norm * param$sigmaX + param$muX
  p.value <- pchisq(Q.Norm1,  df = param$l, ncp=param$d, lower.tail=FALSE, log.p=log.p)
  re <- list(p.value = p.value, Q.Norm=Q.Norm1)  
  return(re)
  
}

Get_Liu_Params_Mod_Lambda <- function(lambda){
  ## Helper function for getting the parameters for the null approximation
  
  c1 <- rep(0,4)
  for(i in 1:4){
    c1[i] <- sum(lambda^i)
  }
  
  re <- Get_Liu_Params_Mod(c1)
  return(re)
}

Get_Liu_PVal.MOD.Lambda.Zero <- function(Q, muQ, muX, sigmaQ, sigmaX, l, d){
  
  
  Q.Norm <- (Q - muQ)/sigmaQ
  Q.Norm1 <- Q.Norm * sigmaX + muX
  
  temp <- c(0.05,10^-10, 10^-20,10^-30,10^-40,10^-50, 10^-60, 10^-70, 10^-80, 10^-90, 10^-100)
  
  out<-qchisq(temp,df = l,ncp=d, lower.tail=FALSE)
  IDX <- max(which(out < Q.Norm1))
  
  pval.msg<-sprintf("Pvalue < %e", temp[IDX])
  return(pval.msg)
}

SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
  
  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  library(SKAT)
  out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0),PACKAGE="SKAT")
  
  out$res <- 1 - out$res
  
  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
  
}

