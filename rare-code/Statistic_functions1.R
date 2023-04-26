LiuParams<-function(a, b, c) {
  re <- list(df = b^2/c, muQ = a, sigmaQ = sqrt(2 * b))
  return(re)
}

pKuonen<-function (x, lambda, delta = rep(0, length(lambda))) 
{ # modified from Thomas Lumley's saddle function in survey package to include noncentrality parameters, which is licensed under GPL-2 | GPL-3
  delta <- delta[lambda > 0]
  lambda <- lambda[lambda > 0]
  if (x <= 0) 
    return(1)
  if (length(lambda) == 1) {
    pchisq(x/lambda, df = 1, ncp = delta, lower.tail = FALSE)
  }
  d <- max(lambda)
  lambda <- lambda/d
  x <- x/d
  k0 <- function(zeta) {
    -sum(log(1 - 2 * zeta * lambda))/2 + sum((delta * lambda * 
                                                zeta)/(1 - 2 * zeta * lambda))
  }
  kprime0 <- function(zeta) {
    sapply(zeta, function(zz) {
      sum(lambda/(1 - 2 * zz * lambda)) + sum((delta * 
                                                 lambda)/(1 - 2 * zz * lambda) + 2 * (delta * 
                                                                                        zz * lambda^2)/(1 - 2 * zz * lambda)^2)
    })
  }
  kpprime0 <- function(zeta) {
    2 * sum(lambda^2/(1 - 2 * zeta * lambda)^2) + sum((4 * 
                                                         delta * lambda^2)/(1 - 2 * zeta * lambda)^2 + 8 * 
                                                        delta * zeta * lambda^3/(1 - 2 * zeta * lambda)^3)
  }
  n <- length(lambda)
  if (any(lambda < 0)) {
    lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
  }
  else if (x > sum(lambda)+sum(delta*lambda)) {
    lmin <- -0.01
  }
  else {
    lmin <- -length(lambda)*max(1+delta)/(2 * x)
  }
  lmax <- min(1/(2 * lambda[lambda > 0])) * 0.99999
  hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, 
                     upper = lmax, tol = 1e-08)$root
  w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
  v <- hatzeta * sqrt(kpprime0(hatzeta))
  if (abs(hatzeta) < 1e-04) 
    NA
  else pnorm(w + log(v/w)/w, lower.tail = FALSE)
}

wuweights <- function(maf){
  ifelse(maf>0, dbeta(maf, 1, 25), 0)
} 

poweights <- function(maf){
 ifelse(maf>0, maf^{-3/4}, 0)
#ifelse(maf>0, 1/sqrt(maf*(1-maf)), 0)
# ifelse(maf>0, 1, 0)
}
### Calculate the covariance matrix of score function for joint effects
GetVar_joint <- function(G, Z, E.inter, mu){
  #G is the genotype for all subjects
  #Z=(1, E, X) is the covariate, including the intercept and the environmental factor
  #E.inter is the normalized environmental factor  g <- MAFs*2
  # mu=E(Y);
  G <- as.matrix(G)
  n <- dim(G)[1]
  m <- dim(G)[2]
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  B11 <- Matrix::bdiag(t(Z)%*%(mu*(1-mu)*Z), n*Sigma/((MAFs*(1-MAFs))%*%t(MAFs*(1-MAFs))))
  B12 <- rbind(cbind(t(Z)%*%(mu*(1-mu))%*%t(g), t(Z)%*%(mu*(1-mu)*E.inter)%*%t(g)), 
               cbind(sum(mu)/(MAFs*(1-MAFs))*Sigma, sum(mu*E.inter)/(MAFs*(1-MAFs))*Sigma) )  
  A <- g%*%t(g)
  B22 <- rbind(cbind(sum(mu)*Sigma, sum(mu*E.inter)*Sigma), cbind(sum(mu*E.inter)*Sigma, sum(E.inter^2*mu)*Sigma))+
    rbind(cbind(sum(mu*(1-mu))*A,  sum(mu*(1-mu)*E.inter)*A), cbind(sum(mu*(1-mu)*E.inter)*A,  sum(mu*(1-mu)*E.inter^2)*A))
  
  eig_B11 <- eigen(B11,symmetric=T)
  id <- which(eig_B11$values>0)
  B11_inv <- eig_B11$vectors[, id] %*% diag(1/eig_B11$values[id])%*%t(eig_B11$vectors[, id])
  
  Sigma.S.star <- B22-t(B12)%*%B11_inv%*%B12
  return(Sigma.S.star) 
}

### Calculate the covariance matrix of score function for interaction effects
GetVar_inter <- function(Y, G, X, nullmodel, weights1){
  
  G <- as.matrix(G)
  E <- X[, 1]
  Z <- cbind(1, X)
  
  n <- length(Y)
  m <- dim(G)[2]
  MAFs <- colMeans(G)/2
  g <- 2*MAFs
  Sigma <- var(G)
  
  beta1 <- as.vector(nullmodel$coefficients$random$id) ## random effects
  beta0 <- as.vector(nullmodel$coefficients$fixed)     ## fixed effects
  G.weight <- t(t(G)*weights1)
  
  #E1 <- colMeans(exp(t(t(G.weight)*beta1)))
  tmp0 <- as.vector(exp(G.weight%*%beta1))
 
  E2 <- colMeans(G*tmp0)
  #E3 <- colMeans(G*exp(t(t(G.weight)*beta1)))
  E4 <- matrix(NA, m, m)
  for (i in 1:m) {
    for(j in i:m){
      if(i==j){
        E4[i, i] <- mean(G[, i]^2*tmp0)
      }else{
        E4[i, j] <- mean(G[, i]*G[, j]*tmp0) 
      }
    }
  }
  E4[lower.tri(E4)] <- t(E4)[lower.tri(E4)]
  
  tmp1 <- as.vector(exp(Z%*%beta0)) 
  tmp2 <- 1+as.vector(exp(Z%*%beta0))*mean(tmp0)
 
  EGY <- tmp1%*%t(E2)/tmp2
  EG <- (sapply(g, rep, n)+tmp1%*%t(E2))/tmp2
  EY <- (tmp2-1)/tmp2
 
  phi <- nullmodel$phi
 # tau <- nullmodel$tau
 # W <- diag(EY*(1-EY))/phi
 # Gw <- W%*%G.weight
 # D <- solve(diag(1/tau, m) + t(G.weight)%*%Gw)
 # Vinv <- W-Gw %*% D %*% t(Gw)
   Vinv <- diag(EY*(1-EY))/phi
 
  a11 <- t(Z)%*%Vinv%*%Z    #t(Z)%*%diag(EY*(1-EY)/phi)%*%Z   
  a12 <- t(Z)%*%(EGY-EG*EY) 
  a13 <- t(Z)%*%(EGY*(1-EY))  
  a14 <- t(Z)%*%(E*(EGY*(1-EY))) 
  
  a23 <- (sum(tmp1/tmp2)*E4-t(EG)%*%EGY)   
  a24 <- (sum(E*tmp1/tmp2)*E4-t(E*EG)%*%EGY) 
  a21 <- t(a12)
  a22 <- sum(1/tmp2)*(Sigma+g%*%t(g))+sum(tmp1/tmp2)*E4-t(EG)%*%EG  
  
  a31 <- t(a13)
  a32 <- t(a23)
  a33 <-  sum(tmp1/tmp2)*E4-t(EGY)%*%EGY                 
  a34 <- sum(E*tmp1/tmp2)*E4-t(E*EGY)%*%EGY   
  
  a41 <- t(a14)
  a42 <- t(a24)
  a43 <- t(a34)
  a44 <- sum(E^2*tmp1/tmp2)*E4-t(E*EGY)%*%(E*EGY)       
  
  B11 <-  rbind(rbind(cbind(a11, a12, a13), cbind(a21, a22, a23)), cbind(a31, a32, a33))
  B12 <- rbind(rbind(a14, a24), a34)
  B22 <- a44   
  eig_B11 <- eigen(B11,symmetric=T)
  id <- which(eig_B11$values>0)
  B11_inv <- eig_B11$vectors[, id] %*% diag(1/eig_B11$values[id])%*%t(eig_B11$vectors[, id])
  Sigma.star <-  B22-t(B12)%*%B11_inv%*%B12  #Sigma.star[1:5, 1:5]
  return(Sigma.star)
}

### Calculate the p values of interaction effect based on aGE method
GetInter_pval <- function(Sigma, U, pow=pow, n.perm=n.perm){
  CovS <- Sigma
  CovS.edecomp <- eigen(CovS, symmetric=T)
  id <- which(CovS.edecomp$values>0)
  Sigma <- t(CovS.edecomp$vectors[, id] %*% (t(CovS.edecomp$vectors[, id]) * CovS.edecomp$values[id])) #V22 is the standard deviation
  
  Ts <- rep(NA, length(pow))
  for (j in 1:length(pow)) {
    if (pow[j] < Inf) {
      Ts[j] <- sum((U)^pow[j])
    }
    else {
      Ts[j] <- max(abs(U))
    }
  }
   
  U0 <- MASS::mvrnorm(n.perm, rep(0, length(U)), Sigma) 
  pPerm0 <- rep(NA,length(pow))
  T0s <- matrix(NA,n.perm,length(pow)) 
  for (j in 1:length(pow)){
    if (pow[j] < Inf){ T0s[,j] = round(rowSums((U0)^pow[j]), digits = 8) }
    if (pow[j] == Inf) {T0s[,j] = round(apply(abs(U0), 1, max), digits = 8) }
  } 
  
  T0s <- T0s[complete.cases(T0s),]
  n.perm <- nrow(T0s)
  P0s <- matrix(NA,n.perm,length(pow))
  for (j in 1:length(pow)){
    pPerm0[j] <- round((sum(abs(Ts[j])<=abs(T0s[1:(n.perm-1),j]))+1)/(n.perm), digits = 8)
    P0s[,j] <- (n.perm-rank(abs(T0s[,j]))+1)/(n.perm)
  }
  
  minp0 <- apply(P0s,1,min)
  fisherp0 <- apply(P0s,1,function(x) sum(log(x)))
  Paspu <- (sum(minp0<=min(pPerm0))+1)/(n.perm+1)
  Pfisher <- (sum(fisherp0<=sum(log(pPerm0)))+1)/(n.perm+1)
  Ts <- c(Ts, min(pPerm0),sum(log(pPerm0)))
  pvs <- c(pPerm0, Paspu,Pfisher)
  names(pvs) <- c(paste("aGEsm", pow, sep = ""), "aGEsm","aGEsm_fisher")
  return(pvs)
}

### Calculate the p values of joint effect based on aGE method
GetJoint_pval <- function(Sigma, U, X, res, G, EG=G.joint, n=n, pow=pow, n.perm=n.perm){
  
  m <- length(U)/2
  U.main <- U[1:m]
  U.int <- U[-c(1:m)]
  Tstar.main <- rep(NA, length(pow))
  Tstar.int <- rep(NA, length(pow))
  Tstar <- rep(NA,length(pow))
  npow.star <- pow
  for (j in 1:length(pow)) {
    if (pow[j] < Inf) {
      Tstar.main[j] <- sum((U.main^pow[j]))
      Tstar.int[j] <- sum((U.int^pow[j]))
    }
    else {
      Tstar.main[j] <- max(abs(U.main))
      Tstar.int[j] <- max(abs(U.int))
      npow.star[j] <- 0
    }
  }
  Tstar=Tstar.main+Tstar.int
  
  ### Use Monte Carlo method to calculate the p values of aGE ###########
  XdE <- X[, -1, drop=F]
  
  CovStar <- Sigma/n
  CovStar.edecomp <- eigen(CovStar,symmetric=T)
  CovStar.ev <- ifelse(zapsmall(CovStar.edecomp$values)>0, CovStar.edecomp$values, 0)
  V11.star <- t(CovStar.edecomp$vectors %*% (t(CovStar.edecomp$vectors) * sqrt(CovStar.ev)))
  
  s <- sample(1:10^5,1)
  set.seed(s)
  pvs <- calculation(method='Simulation', n.perm=n.perm, res=res, XdE=XdE, nZ=m, n=n, Us=Us, V11=V11.star, 
                     pow=pow, Ts_main=Tstar.main,Ts_int=Tstar.int,Ts=Tstar, G=G, EG=G.joint[, -c(1:m)])
  return(pvs)
}

### To test the independence between E and SNPs: 
### Calculate the covariance matrix of score function for genetic main effects
GetInd_pval <- function(Y, G, X, alpha0, n.perm = n.perm){
 
  G <- as.matrix(G)
  X <- as.matrix(X)
  m <- ncol(G)
  
  MAFs <- colMeans(G, na.rm = T)/2
  E.inter <- as.vector(scale(X[, 1], center = T, scale = T))
  id <- which(Y==0)
  cov <- X[id, -1]
  Y <- E.inter[id]
  X <- G[id,]
  ### Test the independence between E and causal SNPs 
  rest <- aSPU::aSPU(Y=Y, X=X, cov = cov, resample = "perm", model = "gaussian", n.perm = n.perm)
  p0 <- rest$pvs["aSPU"]
  #pval <- apply(G[id, ], 2, function(x) kruskal.test(E.inter[id]~x)$p.value) 
  #p0 <- pcauchy(as.numeric(mean(tan((0.5-pval)*pi))), scale = 1, lower.tail = FALSE)
  
  rho <- ifelse(p0 > alpha0, 0, 1)
  re <- list(p0=p0, rho=rho)
  return(re)
}

## Using Monte Carlo method to caculate the p value of rareGE based on retrospective method
GetrareGE.joint_pval <- function(Uz, V, n=n, rho=seq(0, 1, 0.1), n.perm=n.perm){
  
  ### Uz is the score function for joint effect
  ##V=Sigma.z is the covariance matrix of Uz
  m <- length(Uz)/2
  T0 <- p0 <- matrix(NA,nrow = n.perm, ncol=length(rho))
  Q <- p  <- rep(NA, length(rho))
  
  CovStar.edecomp <- eigen(V,symmetric=T)
  CovStar.ev <- ifelse(zapsmall(CovStar.edecomp$values)>0, CovStar.edecomp$values, 0)
  V <- t(CovStar.edecomp$vectors %*% (t(CovStar.edecomp$vectors) * (CovStar.ev)))
  
  U0 <- MASS::mvrnorm(n.perm, rep(0, length(Uz)), V)
  for (j in 1:length(rho)){
    weights <- c(rho[j]*rep(1, m), (1-rho[j])*rep(1, m))
    T0[, j] <- colSums(weights*t(U0^2))  
    Sigma_half <- Get_Matrix_Square.1(V)
    M <- Sigma_half %*% diag(weights) %*% t(Sigma_half)
    out.s <- eigen(M,symmetric=TRUE, only.values = TRUE)
    lambda1 <- out.s$values
    IDX1 <- which(lambda1 >= 0)
    IDX2 <- which(lambda1 > mean(lambda1[IDX1])/100000)
    lambda <- lambda1[IDX2]
    p0[, j] <- Get_PValue.Lambda(lambda,T0[, j])$p.value
    
    Q[j] <- sum(weights*t(Uz^2))
    p[j] <- Get_PValue.Lambda(lambda,Q[j])$p.value
    #Get_Davies_PVal(Q[j], M)$p.value
    
  }
  
  T.min <- apply(p0, 1, min)
  p.min <- (sum(T.min <= min(p))+1)/(n.perm+1)
  re <- list(T.min=T.min, pval=p.min, p0=p)
  return(re)
}

GetrareGE.joint_pval2 <- function(Uz, V, n=n, rho=seq(0, 1, 0.1), n.perm=n.perm){
  
  ### Uz is the score function for joint effect
  ##V=Sigma.z is the covariance matrix of Uz
  
  m <- length(Uz)/2
  Q <- rho*sum(Uz[1:m]^2)+(1-rho)*sum(Uz[-(1:m)]^2)
  lambdas <- NULL
  liuparams <- NULL
  ps <- Ts <- numeric(length(rho))
  V11<- V[1:m, 1:m]   #t(G1)%*%V%*%G1
  V12<- V[1:m, -c(1:m)] #t(G1)%*%V%*%G2
  V22<- V[-c(1:m), -c(1:m)] #t(G2)%*%V%*%G2
  V11V11 <- V11 %*% V11
  V11V12 <- V11 %*% V12
  V12V22 <- V12 %*% V22
  V12V12 <- V12 %*% t(V12)
  V22V22 <- V22 %*% V22
  a1 <- sum(diag(V11))
  a2 <- sum(diag(V22))
  b1 <- sum(diag(V11V11))
  b2 <- sum(diag(V12V12))
  b3 <- sum(diag(V22V22))
  c1 <- sum(diag(V11V11 %*% V11V11))
  c2 <- sum(diag(V11V12 %*% t(V11V12)))
  c3 <- sum(diag(V11V12 %*% t(V12V22)))
  c4 <- sum(diag(V12V12 %*% V12V12))
  c5 <- sum(diag(V12V22 %*% t(V12V22)))
  c6 <- sum(diag(V22V22 %*% V22V22))
  for (i in 1:length(rho)) {
    lam <- rho[i] * a1 + (1-rho[i]) * a2
    lam2 <- rho[i]^2 * b1 + 2 * rho[i] * (1-rho[i]) * b2 + (1-rho[i])^2 * b3
    lam4 <- rho[i]^4 * c1 + 4 * rho[i]^3 * (1-rho[i]) * c2 + rho[i]^2 * (1-rho[i])^2 * (4 * c3 + 2 * c4) + 4 * rho[i] * (1-rho[i])^3 * c5 + (1-rho[i])^4 * c6
    tmp.param <- LiuParams(lam, lam2, lam4)
    ps[i] <- pchisq((Q[i]-tmp.param$muQ)/tmp.param$sigmaQ*sqrt(2*tmp.param$df)+tmp.param$df, tmp.param$df, lower.tail=F)
    liuparams <- c(liuparams, list(tmp.param))
  }
  minp <- min(ps)
  for (i in 1:length(rho)) {
    Ts[i] <- (qchisq(minp, liuparams[[i]]$df, lower.tail=F)-liuparams[[i]]$df)/sqrt(2*liuparams[[i]]$df)*liuparams[[i]]$sigmaQ+liuparams[[i]]$muQ
  }
  
  V.cond <- V22-t(V12)%*%ginv(V11)%*%V12
  eig <- eigen(V.cond, symmetric = TRUE)
  eigval<-eig$values
  D <- diag(eigval)
  diag(D)[zapsmall(diag(D)) > 0] <- 1/sqrt(diag(D)[zapsmall(diag(D)) > 0])
  diag(D)[diag(D) <= 0] <- 0
  meanvec <- D %*% t(eig$vectors) %*% t(V12) %*% ginv(V11)
  
  IDX<-which(rho>=0.999)
  if(length(IDX)>0) {
    rho[IDX]<-0.999
  }
  
  Fcond<-function(x) {
    qtmp<-min((Ts-rho*sum(x^2))/(1-rho))
    ptmp<-pKuonen(qtmp, lambda=eigval, delta=c(meanvec%*%x)^2)
    return(ptmp)
  }
  
  tmpa<-mvrnorm(n.perm, rep(0, nrow(V11)), V11)
  tmpres<-apply(tmpa, 1, Fcond)
  actualp<-mean(tmpres, na.rm=T)
  
  if(length(IDX)>0) {
    rho[IDX]<-1
  }
  
  re <- list(minp=minp, pval=actualp, rho=rho[which.min(ps)], ps=ps)
  return(re)
}

##  caculate the p value of interaction effect based rareGE
GetrareGE.inter_pval <- function(Sigma, U){
  Q <- sum(U^2)
  eig <- eigen(Sigma, symmetric=T)
  evals <- eig$values[which(eig$values>0)]
  pint_ran <- survey::pchisqsum(Q, rep(1,length(evals)), evals, lower.tail=F, method="sad")
  
  #Sigma_half <- Get_Matrix_Square.1(Sigma)
  #M <- Sigma_half %*% t(Sigma_half)
  #pint_ran <- Get_Davies_PVal(Q, M)$p.value
  
  return(pint_ran)
}

Get_joint_HC_glm <- function(data, alpha=0.05){
  
  Y <- data$Y
  X <- data$X
  G <- data$G
  Z <- cbind(1, X)
  
  n <- dim(G)[1]
  m <- dim(G)[2]
  k <- ncol(Z)
  
  E.inter <- as.vector(scale(X[, 1], center = T, scale = T))
  G.joint <- cbind(G, G*E.inter)
  Z.joint <- cbind(1, E.inter)
  m.joint <- dim(G.joint)[2]
  n.joint <- dim(Z.joint)[2]
  
  fit0 <- glm(Y ~ X,family = "binomial")
  mu <- as.vector(fit0$fitted.values)   ## the mean of Y
  res <- Y-mu  
  
  S <- as.vector(t(G.joint)%*%res)
  E <- diag(mu*(1-mu))
  Sigma.S <- t(G.joint)%*%(E-(Z*mu*(1-mu))%*%solve(t(Z)%*%(Z*mu*(1-mu)))%*%t(Z*mu*(1-mu)))%*%G.joint
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  Sigma.S.star <- GetVar_joint(G, Z, E.inter, mu)
  S.star <- as.vector(t(G.joint)%*%Y)-c(g*sum(mu), g*sum(E.inter*mu))
  ##########################################################
  Z.s <- as.vector(abs(S/sqrt(diag(Sigma.S))))
  Sigma.Z.s <- diag(1/sqrt(diag(Sigma.S)))%*%Sigma.S%*%diag(1/sqrt(diag(Sigma.S)))
  p.s <- 2*pnorm(Z.s, mean=0, sd=1, lower.tail = FALSE)
  
  Z.s.star <- as.vector(abs(S.star/sqrt(diag(Sigma.S.star))))
  Sigma.Z.s.star <- diag(1/sqrt(diag(Sigma.S.star)))%*%Sigma.S.star%*%diag(1/sqrt(diag(Sigma.S.star)))
  p.s.star <- 2*pnorm(Z.s.star, mean=0, sd=1, lower.tail = FALSE)
  ###########################################
  rest <- aSPU(X[,1], G, cov = X[,-1], resample = "sim",model = "gaussian", pow = c(1:8, Inf),n.perm = 1000,threshold = 10^5)
  p0 <-rest$pvs["aSPU"]
  pval <- apply(G, 2, function(x) kruskal.test(E.inter~x)$p.value) 
  p.acat <- pcauchy(as.numeric(mean(tan((0.5-pval)*pi))), scale = 1, lower.tail = FALSE)
  
  rho <- ifelse(p0>2*(1-alpha)*alpha, 0, 1)
  rho.roc <- ifelse(p.acat>2*(1-alpha)*alpha, 0, 1)
  rho <- as.numeric(rho)
  rho.roc <- as.numeric(rho.roc)
  #########################################################
  
  re <- list(Z_score=list(Z.s=Z.s, Z.s.star=Z.s.star), 
             p_score=list(p.s=p.s, p.s.star=p.s.star),
             Sigma=list(Sigma.Z.s=Sigma.Z.s, Sigma.Z.s.star=Sigma.Z.s.star),
             rho=rho, rho.roc=rho.roc)
  return(re)
}

Get_main_EB_glmm <- function(data, alpha=0.05, pow = c(1:6, Inf), n.perm =1000, 
                             mainweights=wuweights, interweights=wuweights){
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  n <- dim(G)[1]
  m <- dim(G)[2]
  k <- ncol(X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  EG <- t(t(G)*weights2)*E
  G.weight <- t(t(G)*weights1)
  
  tdat <- data.frame(trait=Y, cbind(E, X[, -1]), G=EG, id=rep(1, n))
  if (k == 1) {
    fixedmod <- "trait ~ E"
  } else {
    fixedmod <- paste("trait ~", paste(names(tdat)[2:(k+1)], collapse = " + "))}
  if (m == 1) {
   randommod <- "~ 0 + G"} else {
     randommod <- paste("~ 0 +", paste(names(tdat)[(k + 2):(1 + k + m)], collapse = " + "))
   }
  
  nullmodel <- MASS::glmmPQL(fixed = as.formula(fixedmod), random = list(id = nlme::pdIdent(as.formula(randommod))), data = tdat, family = "binomial", verbose = F)
  mu <- tryCatch(as.numeric(predict(nullmodel,data=tdat,type='response',level=1)),error=function(c) NA, warning = function(c) NA)
  varc <- tryCatch(nlme::VarCorr(nullmodel), error=function(c) NA, warning = function(c) NA)
  tau <- as.numeric(varc[1, 1])
  phi <- as.numeric(varc[nrow(varc),1])
  res <- Y-mu
  W <- diag(mu*(1-mu))
  Gw <- W%*%EG
  D <- solve(diag(1/tau, m) + t(EG)%*%Gw)
  Vinv <- W-Gw %*% D %*% t(Gw)
  
  #fit0 <- glm(Y ~ X,family = "binomial")
  #mu <- as.vector(fit0$fitted.values)   ## the mean of Y
  #res <- Y-mu      ## Y-mu
  #Vinv <- diag(mu*(1-mu))
  
  ### The genetic main effect score function
  S <- as.vector(t(G.weight)%*%res)
  Sigma.S <- t(G.weight)%*%(Vinv-Vinv%*%Z%*%(solve(t(Z)%*%Vinv%*%Z))%*%t(Z)%*%Vinv)%*%G.weight
  
  S.star <- (as.vector(t(G.weight)%*%Y)-g*sum(mu)*weights1)
  Sigma.star <- diag(weights2)%*%GetVar_main(G, X, mu, Vinv)%*%diag(weights2)
  
  ### aSPU fisher######
  aspu.star <- GetInter_pval(Sigma=Sigma.star, U=S.star, pow=pow, n.perm=n.perm)
  p.raspu <- aspu.star
  
  aspu.s <- GetInter_pval(Sigma=Sigma.S, U=S,pow=pow, n.perm=n.perm)
  p.aspu <- aspu.s 
   
  ### Test the independence between E and causal SNPs 
  re <- GetInd_pval(Y, G, X, alpha, n.perm = n.perm)
  rho <- re$rho; p0 <- re$p0
  
  if(rho==1){
    p.aspu.EB <- p.aspu
  }else{
    t <- (1-p0)*(p.aspu-p.raspu)^2
    k <- t/(t+p.aspu*p.raspu)
    pval <- k*tan((0.5-p.aspu)*pi)+(1-k)*tan((0.5-p.raspu)*pi)
    p.aspu.EB <- pcauchy(as.numeric(pval), scale = 1, lower.tail = FALSE)
  }
  re <- list(aspu=c(p.aspu=p.aspu, p.raspu=p.raspu, p.aspu.EB=p.aspu.EB),rho=rho, p0=p0)
  return(re)
}
################################################################
Get_main_realdata <- function(data, alpha=0.05, pow = c(1:6, Inf), n.perm =1000, 
                              mainweights=wuweights, interweights=wuweights){
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  n <- dim(G)[1]
  m <- dim(G)[2]
  k <- ncol(X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  EG <- t(t(G)*weights2)*E
  G.weight <- t(t(G)*weights1)
  
  tdat <- data.frame(trait=Y, cbind(E, X[, -1]), G=EG, id=rep(1, n))
  if (k == 1) {
    fixedmod <- "trait ~ E"
  } else {
    fixedmod <- paste("trait ~", paste(names(tdat)[2:(k+1)], collapse = " + "))}
  if (m == 1) {
    randommod <- "~ 0 + G"} else {
      randommod <- paste("~ 0 +", paste(names(tdat)[(k + 2):(1 + k + m)], collapse = " + "))
    }
  
  nullmodel <- MASS::glmmPQL(fixed = as.formula(fixedmod), random = list(id = nlme::pdIdent(as.formula(randommod))), data = tdat, family = "binomial", verbose = F)
  mu <- tryCatch(as.numeric(predict(nullmodel,data=tdat,type='response',level=1)),error=function(c) NA, warning = function(c) NA)
  varc <- tryCatch(nlme::VarCorr(nullmodel), error=function(c) NA, warning = function(c) NA)
  tau <- as.numeric(varc[1, 1])
  phi <- as.numeric(varc[nrow(varc),1])
  res <- Y-mu
  W <- diag(mu*(1-mu))
  Gw <- W%*%EG
  D <- solve(diag(1/tau, m) + t(EG)%*%Gw)
  Vinv <- W-Gw %*% D %*% t(Gw)
  
  #fit0 <- glm(Y ~ X,family = "binomial")
  #mu <- as.vector(fit0$fitted.values)   ## the mean of Y
  #res <- Y-mu      ## Y-mu
  #Vinv <- diag(mu*(1-mu))
  
  ### The genetic main effect score function
  S <- as.vector(t(G.weight)%*%res)
  Sigma.S <- t(G.weight)%*%(Vinv-Vinv%*%Z%*%(solve(t(Z)%*%Vinv%*%Z))%*%t(Z)%*%Vinv)%*%G.weight
  
  S.star <- (as.vector(t(G.weight)%*%Y)-g*sum(mu)*weights1)
  Sigma.star <- diag(weights1)%*%GetVar_main(G, X, mu, Vinv)%*%diag(weights1)
  
  ### aSPU fisher######
  aspu.star <- GetInter_pval(Sigma=Sigma.star, U=S.star, pow=pow, n.perm=n.perm)
  p.raspu <- aspu.star; p.raspumin <- min(aspu.star[c("aGEsm1","aGEsm2", "aGEsm_fisher")])
  
  aspu.s <- GetInter_pval(Sigma=Sigma.S, U=S,pow=pow, n.perm=n.perm)
  p.aspu <- aspu.s; p.aspumin <- min(aspu.s[c("aGEsm1","aGEsm2", "aGEsm_fisher")])
  
  n.perm <- 2000
  while(p.raspumin<5/n.perm){
    n.perm <- 10*n.perm
    aspu.star <- GetInter_pval(Sigma=Sigma.star, U=S.star, pow=pow, n.perm=n.perm)
    p.raspu <- aspu.star; p.raspumin <- min(aspu.star[c("aGEsm1","aGEsm2", "aGEsm_fisher")])
    
    if(n.perm==2e6)
      break
  }
  
  n.perm <- 2000
  while(p.aspumin<5/n.perm){
    n.perm <- 10*n.perm
    aspu.s <- GetInter_pval(Sigma=Sigma.S, U=S,pow=pow, n.perm=n.perm)
    p.aspu <- aspu.s; p.aspumin <- min(aspu.s[c("aGEsm1","aGEsm2", "aGEsm_fisher")])
    
    if(n.perm==2e6)
      break
  }
  
  ### Test the independence between E and causal SNPs 
  re <- GetInd_pval(Y, G, X, alpha, n.perm = n.perm)
  rho <- re$rho; p0 <- re$p0
  if(rho==1){
    p.aspu.EB <- p.aspu
  }else{
    t <- (1-p0)*(p.aspu-p.raspu)^2
    k <- t/(t+p.aspu*p.raspu)
    pval <- k*tan((0.5-p.aspu)*pi)+(1-k)*tan((0.5-p.raspu)*pi)
    p.aspu.EB <- pcauchy(as.numeric(pval), scale = 1, lower.tail = FALSE)
  }
  re <- list(aspu=c(p.aspu=p.aspu, p.raspu=p.raspu, p.aspu.EB=p.aspu.EB), rho=rho, p0=p0)
  return(re)
}
######################################################################################
Get_inter_glmm <- function(data, alpha=0.05, pow = c(1:6, Inf), n.perm =1000, 
                              mainweights=mainweights, interweights=interweights){
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  
  n <- dim(G)[1]
  n0 <- length(which(Y==0))
  n1 <- length(which(Y==1))
  m <- dim(G)[2]
  k <- ncol(X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  #E <- X[, 1]
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  EG <- t(t(G)*weights2)*E
  G.weight <- t(t(G)*weights1)
  
  ###  aGE and rareGE ######
  aGE <- aGE::aGE(Y, G.weight, cov=X, model = "binomial", pow = pow, n.perm = n.perm, nonparaE = F)
  p.aGE <- aGE[1:9]
  
  p.rareGE <-  rareGE::INT_RAN(phenotype=Y, genotypes=G, covariates=X, mainweights = mainweights, 
                              interweights = interweights, family = "binomial", binomialimpute = T)
 # p.rareGE <- p.rareGE$pINT_RAN
  #######################################################################################WS
  ### Test the independence between E and causal SNPs 
  re <- GetInd_pval(Y, G, X, alpha, n.perm = n.perm)
  rho <- re$rho; p0 <- re$p0
#######################################################################
  tdat <- data.frame(trait=Y, X, G=G.weight, id=rep(1, n))
  if (k == 1) {
    fixedmod <- "trait ~ E"
  } else {
    fixedmod <- paste("trait ~", paste(names(tdat)[2:(k+1)], collapse = " + "))}
  if (m == 1) {
    randommod <- "~ 0 + G"} else {
      randommod <- paste("~ 0 +", paste(names(tdat)[(k + 2):(1 + k + m)], collapse = " + "))
    }
  
  nullmodel <- MASS::glmmPQL(fixed = as.formula(fixedmod), random = list(id = nlme::pdIdent(as.formula(randommod))), 
                             data = tdat, family = "binomial", verbose = F, control = list(sigma=1))#
  mu <- tryCatch(as.numeric(predict(nullmodel,data=tdat,type='response',level=1)),error=function(c) NA, warning = function(c) NA)
  varc <- tryCatch(nlme::VarCorr(nullmodel), error=function(c) NA, warning = function(c) NA)
  tau <- as.numeric(varc[1, 1])
  phi <- as.numeric(varc[nrow(varc),1])
  nullmodel$phi <- phi
  nullmodel$tau <- tau
  #######################################################
  beta1 <- as.vector(nullmodel$coefficients$random$id) # random effects
  beta0 <- as.vector(nullmodel$coefficients$fixed)     #  fixed effects
  tmp0 <- as.vector(exp(G.weight%*%beta1))
  
  id <- which(is.infinite(tmp0))
  if(length(id)>0){
    tmp0 <- tmp0[-id] 
    G <- G[-id, ]
    G.weight <- G.weight[-id, ]
  }
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  tmp1 <- as.vector(exp(Z%*%beta0)) 
  tmp2 <- 1+as.vector(exp(Z%*%beta0))*mean(tmp0)
  
  EGY <- tmp1%*%t(colMeans(G*tmp0))/tmp2
  S.star <- as.vector(t(EG)%*%Y)-colSums(E*EGY)*weights2
  S.star <- S.star/phi
  Sigma.star <- GetVar_inter(Y, G, X, nullmodel, weights1)
  Sigma.star <- diag(weights2)%*%Sigma.star%*%diag(weights2)
  ######################################################
  ###############IndaGE and IndrareGE ############################# 
  IndaGE <- GetInter_pval(Sigma=Sigma.star, U=S.star, pow=pow, n.perm=n.perm)
  p.IndaGE <- IndaGE 
  names(p.IndaGE) <- c(paste("IndaGEsm", pow, sep = ""), "IndaGEsm","IndaGEsm_fisher") 
  
  p.IndrareGE <- GetrareGE.inter_pval(Sigma.star, S.star)
  
  if(rho==1){
    p.waGE <- p.aGE
    p.wrareGE <- p.rareGE
  }else{
    t1 <- (1-p0)*(p.aGE-p.IndaGE)^2
    k1 <- t1/(t1+p0*p.aGE*p.IndaGE)
    TwaGE <- k1*tan((0.5-p.aGE)*pi)+(1-k1)*tan((0.5-p.IndaGE)*pi)
    p.waGE <- pcauchy(as.numeric(TwaGE), scale = 1, lower.tail = FALSE)
    
    t2 <- (1-p0)*(p.rareGE-p.IndrareGE)^2
    k2 <- t2/(t2+p0*p.rareGE*p.IndrareGE)
    TwrareGE <- k2*tan((0.5-p.rareGE)*pi)+(1-k2)*tan((0.5-p.IndrareGE)*pi)
    p.wrareGE <- pcauchy(as.numeric(TwrareGE), scale = 1, lower.tail = FALSE)
  }
  names(p.waGE) <- c(paste("waGEsm", 1:(length(p.aGE)-2), sep = ""), "waGEsm","waGEsm_fisher") 
  
  re <- list(p.aGE=p.aGE, p.rareGE=p.rareGE, 
             p.IndaGE=p.IndaGE, p.IndrareGE=p.IndrareGE,
             p.waGE=p.waGE, p.wrareGE=p.wrareGE, p0=p0)
  return(re)
}
################################################################
Get_inter_realdata <- function(data, alpha=0.05, pow = c(1:6, Inf), n.perm =1000, 
                              mainweights=wuweights, interweights=wuweights){
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
   
  n <- dim(G)[1]
  m <- dim(G)[2]
  k <- ncol(X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  EG <- t(t(G)*weights2)*E
  G.weight <- t(t(G)*weights1)
  ###  aGE and rareGE ######
  aGE <- aGE::aGE(Y, G.weight, cov=X, model = "binomial", pow = pow, n.perm = n.perm, nonparaE = F)
  p.aGE <- aGE[1:9]; p.aGEfisher <- p.aGE[9]
  
  p.rareGE <-  rareGE::INT_RAN(phenotype=Y, genotypes=G, covariates=X, mainweights = mainweights, 
                               interweights = interweights, family = "binomial", binomialimpute = T)
  
  n.perm <- 1000
  while(p.aGEfisher<2.8/n.perm){
    n.perm <- 10*n.perm
    aGE <- aGE::aGE(Y, G.weight, cov=X, model = "binomial", pow = pow, n.perm = n.perm, nonparaE = F)
    p.aGE <- aGE[1:9]; p.aGEfisher <- p.aGE[9]
    
    if(n.perm==1e6)
      break
  }
  ### Test the independence between E and causal SNPs 
  re <- GetInd_pval(Y, G, X, alpha, n.perm = n.perm)
  rho <- re$rho; p0 <- re$p0
  
  tdat <- data.frame(trait=Y, X, G=G.weight, id=rep(1, n))
  if (k == 1) {
    fixedmod <- "trait ~ E"
  } else {
    fixedmod <- paste("trait ~", paste(names(tdat)[2:(k+1)], collapse = " + "))}
  if (m == 1) {
    randommod <- "~ 0 + G"} else {
      randommod <- paste("~ 0 +", paste(names(tdat)[(k + 2):(1 + k + m)], collapse = " + "))
    }
  
  nullmodel <- MASS::glmmPQL(fixed = as.formula(fixedmod), random = list(id = nlme::pdIdent(as.formula(randommod))), 
                             data = tdat, family = "binomial", verbose = F, control = list(sigma=1))#
  mu <- tryCatch(as.numeric(predict(nullmodel,data=tdat,type='response',level=1)),error=function(c) NA, warning = function(c) NA)
  varc <- tryCatch(nlme::VarCorr(nullmodel), error=function(c) NA, warning = function(c) NA)
  tau <- as.numeric(varc[1, 1])
  phi <- as.numeric(varc[nrow(varc),1])
  nullmodel$phi <- phi
  nullmodel$tau <- tau
  #######################################################
  beta1 <- as.vector(nullmodel$coefficients$random$id) # random effects
  beta0 <- as.vector(nullmodel$coefficients$fixed)     #  fixed effects
  tmp0 <- as.vector(exp(G.weight%*%beta1))
  
  id <- which(is.infinite(tmp0))
  if(length(id)>0){
    tmp0 <- tmp0[-id] 
    G <- G[-id, ]
    G.weight <- G.weight[-id, ]
  }
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  tmp1 <- as.vector(exp(Z%*%beta0)) 
  tmp2 <- 1+as.vector(exp(Z%*%beta0))*mean(tmp0)
  
  EGY <- tmp1%*%t(colMeans(G*tmp0))/tmp2
  S.star <- as.vector(t(EG)%*%Y)-colSums(E*EGY)*weights2
  S.star <- S.star/phi
  Sigma.star <- GetVar_inter(Y, G, X, nullmodel, weights1)
  Sigma.star <- diag(weights2)%*%Sigma.star%*%diag(weights2)
  #############################################################
  ###############IndaGE and IndrareGE ############################# 
  IndaGE <- GetInter_pval(Sigma=Sigma.star, U=S.star, pow=pow, n.perm=n.perm)
  p.IndaGE <- IndaGE; p.IndaGEfisher <- IndaGE[9] 
  names(p.IndaGE) <- c(paste("IndaGEsm", pow, sep = ""), "IndaGEsm","IndaGEsm_fisher") 
  
  p.IndrareGE <- GetrareGE.inter_pval(Sigma.star, S.star)
  
  n.perm <- 1000
  while(p.IndaGEfisher<2.8/n.perm){
    n.perm <- 10*n.perm
    p.IndaGE <- GetInter_pval(Sigma=Sigma.star, U=S.star, pow=pow, n.perm=n.perm)
    p.IndaGEfisher <- p.IndaGE[9]
    
    if(n.perm==1e6)
      break
  }
 
  if(rho==1){
    p.waGE <- p.aGE
    p.wrareGE <- p.rareGE
  }else{
    t1 <- (1-p0)*(p.aGE-p.IndaGE)^2
    k1 <- t1/(t1+p0*p.aGE*p.IndaGE)
    TwaGE <- k1*tan((0.5-p.aGE)*pi)+(1-k1)*tan((0.5-p.IndaGE)*pi)
    p.waGE <- pcauchy(as.numeric(TwaGE), scale = 1, lower.tail = FALSE)
    
    t2 <- (1-p0)*(p.rareGE-p.IndrareGE)^2
    k2 <- t2/(t2+p0*p.rareGE*p.IndrareGE)
    TwrareGE <- k2*tan((0.5-p.rareGE)*pi)+(1-k2)*tan((0.5-p.IndrareGE)*pi)
    p.wrareGE <- pcauchy(as.numeric(TwrareGE), scale = 1, lower.tail = FALSE)
  }
  
  names(p.waGE) <- c(paste("waGEsm", 1:(length(p.aGE)-2), sep = ""), "waGEsm","waGEsm_fisher") 
  
  re <- list(p.aGE=p.aGE, p.rareGE=p.rareGE, 
             p.IndaGE=p.IndaGE, p.IndrareGE=p.IndrareGE,
             p.waGE=p.waGE, p.wrareGE=p.wrareGE, p0=p0)
  return(re)
}

#########################################################
Get_joint_glm <- function(data, alpha=0.05, pow = c(1:6, Inf), n.perm =1000, 
                             mainweights=mainweights, interweights=interweights){
  Y <- data$Y
  X <- data$X
  G <- as.matrix(data$G)
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  n <- dim(G)[1]
  m <- dim(G)[2]
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  G.weight <- t(t(G)*weights1)
  
  weights <- c(weights1, weights2)
  G.joint <- cbind(G, G*E)
  G.joint.weight <- t(t(G.joint)*weights)
  
  fit0 <- glm(Y ~ X,family = "binomial")
  mu <- as.vector(fit0$fitted.values)   ## the mean of Y
  res <- Y-mu  
  Sigma.star <- diag(weights)%*%GetVar_joint(G, Z, E, mu)%*%diag(weights)
  S.star <- as.vector(t(G.joint.weight)%*%Y)-c(g*sum(mu), g*sum(E*mu))*weights
  
  ################# IndaGE and IndrareGE################################
  IndaGE <- GetJoint_pval(Sigma=Sigma.star, U=S.star, X=X, res=res, G=G, EG=G.joint.weight, n=n, pow=pow, n.perm=n.perm)
  p.IndaGE <- IndaGE$pvs
  names(p.IndaGE) <- c(paste("IndaGEsm", pow, sep = ""), "IndaGEsm","IndaGEsm_fisher") 
  
  ### aGE and rareGE ######
  #IndrareGE <- GetrareGE.joint_pval(Uz=S.star, V=Sigma.star, n=n, rho=seq(0, 1, 0.1), n.perm=n.perm)
 # p.IndrareGE <- IndrareGE$pval
  
  IndrareGE <- GetrareGE.joint_pval2(Uz=S.star, V=Sigma.star, n=n, rho=seq(0, 1, 0.1), n.perm=n.perm)
  p.IndrareGE <- IndrareGE$pval
  
  p.aGE <- aGE::aGE.joint(Y, G.weight, cov=X, model="binomial", pow=pow, n.perm=n.perm, method="Simulation")
 
  rareGE <-  rareGE::JOINT(phenotype=Y, genotypes=G, covariates=X, mainweights = mainweights, 
                            interweights = interweights, family = "binomial", binomialimpute = T, 
                            rho = seq(0, 1, by = 0.1), B = n.perm)
  p.rareGE <- rareGE$pJOINT
################################################
  re <- GetInd_pval(Y, G, X, alpha, n.perm = n.perm)
  rho <- re$rho; p0 <- re$p0
  if(rho==1){
    p.waGE <- p.aGE
    p.wrareGE <- p.rareGE
  }else{
    t1 <- (1-p0)*(p.aGE-p.IndaGE)^2
    k1 <- t1/(t1+p0*p.aGE*p.IndaGE)
    TwaGE <- k1*tan((0.5-p.aGE)*pi)+(1-k1)*tan((0.5-p.IndaGE)*pi)
    p.waGE <- pcauchy(as.numeric(TwaGE), scale = 1, lower.tail = FALSE)
    
    t2 <- (1-p0)*(p.rareGE-p.IndrareGE)^2
    k2 <- t2/(t2+p0*p.rareGE*p.IndrareGE)
    TwrareGE <- k2*tan((0.5-p.rareGE)*pi)+(1-k2)*tan((0.5-p.IndrareGE)*pi)
    p.wrareGE <- pcauchy(as.numeric(TwrareGE), scale = 1, lower.tail = FALSE)
  }
  re <- list(p.aGE=p.aGE, p.rareGE=p.rareGE, 
             p.IndaGE=p.IndaGE, p.IndrareGE=p.IndrareGE,
             p.waGE=p.waGE, p.wrareGE=p.wrareGE, p0=p0)
  return(re)
}

Get_joint_realdata <- function(data, alpha=0.05, pow = c(1:6, Inf), n.perm =1000, 
                               mainweights=wuweights, interweights=wuweights){
  Y <- data$Y
  X <- data$X
  G <- as.matrix(data$G)
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  n <- dim(G)[1]
  m <- dim(G)[2]
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  G.weight <- t(t(G)*weights1)
  
  weights <- c(weights1, weights2)
  G.joint <- cbind(G, G*E)
  G.joint.weight <- t(t(G.joint)*weights)
  
  fit0 <- glm(Y ~ X,family = "binomial")
  mu <- as.vector(fit0$fitted.values)   ## the mean of Y
  res <- Y-mu  
  Sigma.star <- diag(weights)%*%GetVar_joint(G, Z, E, mu)%*%diag(weights)
  S.star <- as.vector(t(G.joint.weight)%*%Y)-c(g*sum(mu), g*sum(E*mu))*weights
  
  ################# IndaGE and IndrareGE################################
  IndaGE <- GetJoint_pval(Sigma=Sigma.star, U=S.star, X=X, res=res, G=G, EG=G.joint.weight, n=n, pow=pow, n.perm=n.perm)
  p.IndaGE <- IndaGE$pvs; p.IndaGEfisher <- p.IndaGE[10]
  names(p.IndaGE) <- c(paste("IndaGEsm", pow, sep = ""), "IndaGEsm","IndaGEsm_fisher") 
  
  IndrareGE <- GetrareGE.joint_pval2(Uz=S.star, V=Sigma.star, n=n, rho=seq(0, 1, 0.1), n.perm=n.perm)
  p.IndrareGE <- IndrareGE$pval
 ################aGE and rareGE ############################ 
  p.aGE <- aGE::aGE.joint(Y, G.weight, cov=X, model="binomial", pow=pow, n.perm=n.perm, method="Simulation")
  p.aGEfisher <- p.aGE[10]
  
  rareGE <-  rareGE::JOINT(phenotype=Y, genotypes=G, covariates=X, mainweights = mainweights, 
                           interweights = interweights, family = "binomial", binomialimpute = T, 
                           rho = seq(0, 1, by = 0.1), B = n.perm)
  p.rareGE <- rareGE$pJOINT
  #############################################################
  n.perm <- 1000
  while(p.aGEfisher <2.8/n.perm){
    n.perm <- 10*n.perm
    p.aGE <- aGE::aGE.joint(Y, G.weight, cov=X, model="binomial", pow=pow, n.perm=n.perm, method="Simulation")
    p.aGEfisher <- p.aGE[10]
    
    if(n.perm==1e6)
      break
  }
  
  n.perm <- 1000
  while(p.IndaGEfisher < 2.8/n.perm){
    n.perm <- 10*n.perm
    IndaGE <- GetJoint_pval(Sigma=Sigma.star, U=S.star, X=X, res=res, G=G, EG=G.joint.weight, n=n, pow=pow, n.perm=n.perm)
    p.IndaGE <- IndaGE$pvs; p.IndaGEfisher <- p.IndaGE[10]
    
    if(n.perm==1e6)
      break
  }
  
  n.perm <- 1000
  while(p.rareGE<2.8/n.perm){
    n.perm <- 10*n.perm
    rareGE <-  rareGE::JOINT(phenotype=Y, genotypes=G, covariates=X, mainweights = mainweights, 
                             interweights = interweights, family = "binomial", binomialimpute = T, 
                             rho = seq(0, 1, by = 0.1), B = n.perm)
    p.rareGE <- rareGE$pJOINT
    
    if(n.perm==1e6)
      break
  }
  
  n.perm <- 1000
  while(p.IndrareGE<2.8/n.perm){
    n.perm <- 10*n.perm
    IndrareGE <- GetrareGE.joint_pval2(Uz=S.star, V=Sigma.star, n=n, rho=seq(0, 1, 0.1), n.perm=n.perm)
    p.IndrareGE <- IndrareGE$pval
    
    if(n.perm==1e6)
      break
  }
  
  re <- GetInd_pval(Y, G, X, alpha, n.perm = n.perm)
  rho <- re$rho; p0 <- re$p0
  if(rho==1){
    p.waGE <- p.aGE
    p.wrareGE <- p.rareGE
  }else{
    t1 <- (1-p0)*(p.aGE-p.IndaGE)^2
    k1 <- t1/(t1+p0*p.aGE*p.IndaGE)
    TwaGE <- k1*tan((0.5-p.aGE)*pi)+(1-k1)*tan((0.5-p.IndaGE)*pi)
    p.waGE <- pcauchy(as.numeric(TwaGE), scale = 1, lower.tail = FALSE)
    
    t2 <- (1-p0)*(p.rareGE-p.IndrareGE)^2
    k2 <- t2/(t2+p0*p.rareGE*p.IndrareGE)
    TwrareGE <- k2*tan((0.5-p.rareGE)*pi)+(1-k2)*tan((0.5-p.IndrareGE)*pi)
    p.wrareGE <- pcauchy(as.numeric(TwrareGE), scale = 1, lower.tail = FALSE)
  }
  re <- list(p.aGE=p.aGE, p.rareGE=p.rareGE, 
             p.IndaGE=p.IndaGE, p.IndrareGE=p.IndrareGE,
             p.waGE=p.waGE, p.wrareGE=p.wrareGE, p0=p0)
  return(re)
}


 