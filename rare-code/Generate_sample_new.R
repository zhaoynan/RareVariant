args=commandArgs(trailingOnly=TRUE)
M = as.numeric(args[1])
n0 = as.numeric(args[2]) # of controls/cases
nSNP = as.numeric(args[3]) # of non-casual SNPs
rho = as.numeric(args[4]) #AR1(rho) between casual SNPs and non-causal SNPs
eta = as.numeric(args[5]) # the degree of correlation between E and causal SNPs
test = args[6] ## null or alter
set = args[7]  ## null--joint/inter  
variant = args[8]
######################## Generate data ##########################################
set.seed(1234)
library(Matrix)

#variant <- "RV"
dir0 <- paste0("~/Rare-variant", "/", "Inter_test", "/", variant,"/", test)
setwd(dir0)

############## parameters settup #######
nSNP0 <- 8
nSNP <- 50
ORg <- ORgx <-  rep(1, nSNP0)
if(test=="null"){
  
  if(set!="joint"){
    ORg <- runif(nSNP0, 0, 2)
  }
  
}else{
  if(variant=="RV"){
    # interaction test
    #ORgx[id.causal] <- runif(nSNP0, 0.4, 1.5)  
    #ORg[id.causal] <- runif(nSNP0, 0, 2)
    # joint test
    #ORgx[id.causal] <- runif(nSNP0, 0.5, 1.3) 
    #ORg[id.causal] <- runif(nSNP0, 0, 2)
    # compare joint test and interaction test
    ORgx <- runif(nSNP0, 0.6, 1.2) #runif(nSNP0, 0.25, 1.5)
    ORg <- rep(1, nSNP0)
  }else{
    # interaction test
    #ORgx[id.causal] <- runif(nSNP0, 0.7, 1.3) 
    #ORg[id.causal] <- runif(nSNP0, 0, 2)
    # joint test  
    #ORgx[id.causal] <- runif(nSNP0, 0.65, 1.3)
    #ORg[id.causal] <- runif(nSNP0, 0.9, 1.1)
    # compare joint test and interaction test
    ORgx <- c(0.7, 1.1, 1.1, 1.1, 1.1, 1.1, 0.7, 0.7)
    ORg <- rep(1, nSNP0)
  }
}
ORx <- c(0.25, 0.05, 0.5, 0.1)
#################### sample size and prevalence and other parameters settup #########
n1 <- n0
f0 <- 0.01
##################################################
if(variant == "CV-RV"){
  MAF0slow=0.001; MAF0sup=0.45
  MAFslow=0.001; MAFsup=0.45
}else{
  MAF0slow=0.001; MAF0sup=0.05
  MAFslow=0.001; MAFsup=0.05
}

MAF0s <- runif(nSNP0, MAF0slow, MAF0sup)
MAFs <- runif(nSNP, MAFslow, MAFsup)

R0 <- matrix(1, nrow=nSNP0, ncol=nSNP0)
for(i in 1:nSNP0)
  for(j in 1:nSNP0)
    if (i!=j) R0[i, j] <- rho^(abs(i-j))

R0noise <- matrix(1, nrow=nSNP, ncol=nSNP)
for(i in 1:nSNP)
  for(j in 1:nSNP)
    if (i!=j) R0noise[i, j] <- rho^(abs(i-j))
#################################################
para <- list(ORx=ORx, ORg=ORg, ORgx=ORgx, eta=eta, n0=n0, n1=n1, nSNP=nSNP, MAF0s=MAF0s, MAFs=MAFs, R0=R0, R0noise=R0noise, f0)
if(set=="inter"){
  simAR1RareSNP.inter <- function(para){
    ORg <- para$ORg
    ORx <- para$ORx
    ORgx <- para$ORgx
    
    n0 <- para$n0
    n1 <- para$n1
    
    nSNP <- para$nSNP
    R0 <- para$R0
    R0noise <- para$R0noise
    ####################################
    nSNP0 <- length(ORg)  # The number of csRV
    m <- nSNP+nSNP0  # the number of all SNPs
    n <- n1+n0
    
    betag <- log(ORg)
    betax <- ORx
    betagx <- log(ORgx)
    beta <- c(betax, betag, betagx)
    alpha0 <- log(f0/(1-f0))
    
    X1 <- rnorm(n, mean = 30, sd=5.7)
    X2 <- rbinom(n, 1, prob = 0.5)
    X3 <- rnorm(n, mean = 160, sd=9.4)
    
    Y <- rep(0, n); Y[(n0+1):n] <- 1
    # sample controls:
    eigen.R <- eigen(R0, symmetric = T)
    R1 <- eigen.R$vectors%*% diag(sqrt(eigen.R$values))
    MAF0s <- para$MAF0s 
    cutoff0 <- qnorm(MAF0s)
    
    G0 <- matrix(0, nrow = n, ncol = nSNP0)
    X <- rep(NA, n)
    
    i <- 1
    while (i <= n0) {
      r0 <- rnorm(nSNP0, 0, 1)
      r1 <- R1%*%r0 
      r2 <- ifelse(r1<cutoff0, 1, 0)
      
      r0 <- rnorm(nSNP0, 0, 1)
      r1 <- R1%*%r0 
      r3 <- ifelse(r1<cutoff0, 1, 0)
      
      r4 <- r2+r3 ## RV
      
      E <- sum(r4)*log(eta)+rnorm(1, 0, 1)
      ge <-  E*r4 
      
      tmp <- exp(alpha0+sum(betagx*ge))
      pr <- tmp/(1+tmp)
     
      Y1 <- sample(c(0, 1), 1, prob = c(1-pr, pr))
      if(Y1==0){
        G0[i, ] <- r4
        X[i] <- E
        i <- i+1
      }
    }
    
    #sampling cases:
    while(i <= n){
      r0 <- rnorm(nSNP0, 0, 1)
      r1 <- R1%*%r0 
      r2 <- ifelse(r1<cutoff0, 1, 0)
      
      r0 <- rnorm(nSNP0, 0, 1)
      r1 <- R1%*%r0 
      r3 <- ifelse(r1<cutoff0, 1, 0)
      
      r4 <- r2+r3 ## RV
      
      E <- sum(r4)*log(eta)+rnorm(1, 0, 1)
      ge <- E*r4
      tmp <- exp(alpha0+sum(betagx*ge))#exp(alpha0+sum(c(betag, betagx)*c(r4, ge)))
      pr <- tmp/(1+tmp)
      Y1 <- sample(c(0, 1), 1, prob = c(1-pr, pr))
      if(Y1==1){
        G0[i, ] <- r4
        X[i] <- E
        i <- i+1
      }
    }
    if (nSNP>0){
      ## generate neutral SNPs:
      eigen.Rnoise <- eigen(R0noise, symmetric = T)
      R1noise <- eigen.Rnoise$vectors%*% diag(sqrt(eigen.Rnoise$values))
      MAFs <- para$MAFs 
      cutoff <- qnorm(MAFs)
      
      Gnoise <- matrix(0, nrow=n, ncol = nSNP)
      for(i in 1:n){
        r0 <- rnorm(nSNP, 0, 1) #: X0 ~ MVN(0, I)
        r1 <- R1noise %*% r0   #: X1 ~ MVN(0, R)
        r2 <- ifelse(r1 < cutoff, 1, 0)
        
        r0 <- rnorm(nSNP, 0, 1) #: X0 ~ MVN(0, I)
        r1 <- R1noise %*% r0   #: X1 ~ MVN(0, R)
        r3 <- ifelse(r1<cutoff, 1, 0)
        
        r4 <- r2+ r3
        Gnoise[i, ] <- r4
      }
      
      Gall <- cbind(G0, Gnoise)
    } else Gall <- G0
    
    re <- list(Y=Y, G=Gall, X=cbind(X, X1, X2, X3))
    return(re) 
  }
  
  fun <- function(i){
    Data <- simAR1RareSNP.inter(para)
    return(Data)
  }
  
  library(parallel)
  system.time(data <- mclapply(1:M, fun, mc.cores=35))
  data$OR <- c(ORg, ORx, ORgx)
  save(data, file=paste0(dir0, "/", n0, test, "_", set, "_", nSNP, "_", rho, "_", eta, ".RData"))
  rm(list = ls())
  
}else{
  simAR1RareSNP.joint <- function(para){
    ORg <- para$ORg
    ORx <- para$ORx
    ORgx <- para$ORgx
    
    n0 <- para$n0
    n1 <- para$n1
    
    nSNP <- para$nSNP
    R0 <- para$R0
    R0noise <- para$R0noise
    ####################################
    nSNP0 <- length(ORg)  # The number of csRV
    m <- nSNP+nSNP0  # the number of all SNPs
    n <- n1+n0
    
    betag <- log(ORg)
    betax <- ORx
    betagx <- log(ORgx)
    beta <- c(betax, betag, betagx)
    alpha0 <- log(f0/(1-f0))
    
    X1 <- rnorm(n, mean = 30, sd=5.7)
    X2 <- rbinom(n, 1, prob = 0.5)
    X3 <- rnorm(n, mean = 160, sd=9.4)
    
    Y <- rep(0, n); Y[(n0+1):n] <- 1
    # sample controls:
    eigen.R <- eigen(R0, symmetric = T)
    R1 <- eigen.R$vectors%*% diag(sqrt(eigen.R$values))
    MAF0s <- para$MAF0s 
    cutoff0 <- qnorm(MAF0s)
    
    G0 <- matrix(0, nrow = n, ncol = nSNP0)
    X <- rep(NA, n)
    
    i <- 1
    while (i <= n0) {
      r0 <- rnorm(nSNP0, 0, 1)
      r1 <- R1%*%r0 
      r2 <- ifelse(r1<cutoff0, 1, 0)
      
      r0 <- rnorm(nSNP0, 0, 1)
      r1 <- R1%*%r0 
      r3 <- ifelse(r1<cutoff0, 1, 0)
      
      r4 <- r2+r3 ## RV
      
      E <- sum(r4)*log(eta)+rnorm(1, 0, 1)
      ge <-  E*r4 
      
      tmp <- exp(alpha0+sum(c(betag, betagx)*c(r4, ge)))
      pr <- tmp/(1+tmp)
      
      Y1 <- sample(c(0, 1), 1, prob = c(1-pr, pr))
      if(Y1==0){
        G0[i, ] <- r4
        X[i] <- E
        i <- i+1
      }
    }
    
    #sampling cases:
    while(i <= n){
      r0 <- rnorm(nSNP0, 0, 1)
      r1 <- R1%*%r0 
      r2 <- ifelse(r1<cutoff0, 1, 0)
      
      r0 <- rnorm(nSNP0, 0, 1)
      r1 <- R1%*%r0 
      r3 <- ifelse(r1<cutoff0, 1, 0)
      
      r4 <- r2+r3 ## RV
      
      E <- sum(r4)*log(eta)+rnorm(1, 0, 1)
      ge <- E*r4
      tmp <- exp(alpha0+sum(c(betag, betagx)*c(r4, ge)))
      pr <- tmp/(1+tmp)
      Y1 <- sample(c(0, 1), 1, prob = c(1-pr, pr))
      if(Y1==1){
        G0[i, ] <- r4
        X[i] <- E
        i <- i+1
      }
    }
    if (nSNP>0){
      ## generate neutral SNPs:
      eigen.Rnoise <- eigen(R0noise, symmetric = T)
      R1noise <- eigen.Rnoise$vectors%*% diag(sqrt(eigen.Rnoise$values))
      MAFs <- para$MAFs 
      cutoff <- qnorm(MAFs)
      
      Gnoise <- matrix(0, nrow=n, ncol = nSNP)
      for(i in 1:n){
        r0 <- rnorm(nSNP, 0, 1) #: X0 ~ MVN(0, I)
        r1 <- R1noise %*% r0   #: X1 ~ MVN(0, R)
        r2 <- ifelse(r1 < cutoff, 1, 0)
        
        r0 <- rnorm(nSNP, 0, 1) #: X0 ~ MVN(0, I)
        r1 <- R1noise %*% r0   #: X1 ~ MVN(0, R)
        r3 <- ifelse(r1<cutoff, 1, 0)
        
        r4 <- r2+ r3
        Gnoise[i, ] <- r4
      }
      
      Gall <- cbind(G0, Gnoise)
    } else Gall <- G0
    
    re <- list(Y=Y, G=Gall, X=cbind(X, X1, X2, X3))
    return(re) 
  }
  
  fun <- function(i){
    Data <- simAR1RareSNP.joint(para)
    return(Data)
  }
  
  library(parallel)
  system.time(data <- mclapply(1:M, fun, mc.cores=35))
  data$OR <- c(ORg, ORx, ORgx)
  save(data, file=paste0(dir0, "/", n0, test, "_", set, "_", nSNP, "_", rho, "_", eta, ".RData"))
  rm(list = ls())
  
}
