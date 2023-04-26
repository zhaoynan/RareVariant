args=commandArgs(trailingOnly=TRUE)
M = as.numeric(args[1]) # number of replications
n0 = as.numeric(args[2])
eta = as.numeric(args[3])
test = args[4]
set = args[5]
variant = args[6]
  
dir0 <- paste0("~/Rare-variant", "/", "Inter_test", "/", variant,"/", test)
setwd(dir0)

candi_SNPs <- as.vector(read.table(file = "~/PretermData/PLINK_filtered_dataset/SEC16B.txt", header = F))$V1
library(LDlinkR)
# caculate the LD matrix based on 1000 Genomes Project population
LDinfo <- LDmatrix(snps = candi_SNPs,  pop = "EUR", r2d = "r2",  token = '35deec53ae3c',file =FALSE)

colnames(LDinfo) <- NULL
LD <- as.matrix(LDinfo[, -1])
#set causal snps
nSNP <- dim(LD)[1]
# Generate genotype data
library(Matrix)
nSNP0 <- 8 
id.causal <- 1:nSNP0

ORg <- ORgx <- rep(1, nSNP)

set.seed(1234)
if(test=="null"){
  
  if(set!="joint"){
    ORg[id.causal] <- runif(nSNP0, 0, 2)
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
    ORgx[id.causal] <- runif(nSNP0, 0.6, 1.2) #runif(nSNP0, 0.25, 1.5)
  }else{
    # interaction test
    #ORgx[id.causal] <- runif(nSNP0, 0.7, 1.3) 
    #ORg[id.causal] <- runif(nSNP0, 0, 2)
    # joint test  
    #ORgx[id.causal] <- runif(nSNP0, 0.65, 1.3)
    #ORg[id.causal] <- runif(nSNP0, 0.9, 1.1)
    # compare joint test and interaction test
    ORgx[id.causal] <- c(0.7, 1.1, 1.1, 1.1, 1.1, 1.1, 0.7, 0.7)
  }
}

ORx <- c(0.25, 0.05, 0.5, 0.1)

n1 <- n0
n <- n0+n1
f0 <- 0.01

####################################
if(variant == "CV-RV"){
  MAF0slow=0.001; MAF0sup=0.45
}else{
  MAF0slow=0.001; MAF0sup=0.05
}
MAF0s <- runif(nSNP, MAF0slow, MAF0sup)

para <- list(ORx=ORx, ORg=ORg, ORgx=ORgx, eta=eta, n0=n0, n1=n1, nSNP=nSNP, MAF0s=MAF0s, R=LD, f0)
if(set=="inter"){
  simLDRareSNP.inter <- function(para){
    ORg <- para$ORg
    ORx <- para$ORx
    ORgx <- para$ORgx
    
    n0 <- para$n0
    n1 <- para$n1
    
    nSNP <- para$nSNP
    R0 <- para$R 
    ####################################
    n <- n1+n0
    
    betag <- log(ORg)
    betax <- ORx
    betagx <- log(ORgx)
    nSNP0 <- length(which(betag!=0))
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
    
    G0 <- matrix(0, nrow = n, ncol = nSNP)
    X <- rep(NA, n)
    
    i <- 1
    while (i <= n0) {
      r0 <- rnorm(nSNP, 0, 1)
      r1 <- R1%*%r0 
      r2 <- ifelse(r1<cutoff0, 1, 0)
      
      r0 <- rnorm(nSNP, 0, 1)
      r1 <- R1%*%r0 
      r3 <- ifelse(r1<cutoff0, 1, 0)
      
      r4 <- r2+r3 ## RV
      
      if(nSNP0==0){
        E <- sum(r4)*log(eta)+rnorm(1, 0, 1)
      }else{
        E <- sum(r4[1:nSNP0])*log(eta)+rnorm(1, 0, 1)
      }
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
      r0 <- rnorm(nSNP, 0, 1)
      r1 <- R1%*%r0 
      r2 <- ifelse(r1<cutoff0, 1, 0)
      
      r0 <- rnorm(nSNP, 0, 1)
      r1 <- R1%*%r0 
      r3 <- ifelse(r1<cutoff0, 1, 0)
      
      r4 <- r2+r3 ## RV
      
      if(nSNP0==0){
        E <- sum(r4)*log(eta)+rnorm(1, 0, 1)
      }else{
        E <- sum(r4[1:nSNP0])*log(eta)+rnorm(1, 0, 1)
      }
      ge <- E*r4
      tmp <- exp(alpha0+sum(betagx*ge))
      pr <- tmp/(1+tmp)
      Y1 <- sample(c(0, 1), 1, prob = c(1-pr, pr))
      if(Y1==1){
        G0[i, ] <- r4
        X[i] <- E
        i <- i+1
      }
    }
 
    re <- list(Y=Y, G=G0, X=cbind(X, X1, X2, X3))
    return(re) 
  }
  
  fun <- function(i){
    Data <- simLDRareSNP.inter(para)
    return(Data)
  }
  
  library(parallel)
  system.time(data <- mclapply(1:M, fun, mc.cores=40))
  save(data, file=paste0(dir0, "/", n0, test, "_", set, "_", nSNP, "_", "LD", "_", eta, ".RData"))
  rm(list = ls())
}else{
  
  simLDRareSNP.joint <- function(para){
    ORg <- para$ORg
    ORx <- para$ORx
    ORgx <- para$ORgx
    
    n0 <- para$n0
    n1 <- para$n1
    
    nSNP <- para$nSNP
    R0 <- para$R 
    ####################################
    n <- n1+n0
    
    betag <- log(ORg)
    betax <- ORx
    betagx <- log(ORgx)
    alpha0 <- log(f0/(1-f0))
    nSNP0 <- length(which(betag!=0))
      
    X1 <- rnorm(n, mean = 30, sd=5.7)
    X2 <- rbinom(n, 1, prob = 0.5)
    X3 <- rnorm(n, mean = 160, sd=9.4)
    
    Y <- rep(0, n); Y[(n0+1):n] <- 1
    # sample controls:
    eigen.R <- eigen(R0, symmetric = T)
    R1 <- eigen.R$vectors%*% diag(sqrt(eigen.R$values))
    MAF0s <- para$MAF0s 
    cutoff0 <- qnorm(MAF0s)
    
    G0 <- matrix(0, nrow = n, ncol = nSNP)
    X <- rep(NA, n)
    
    i <- 1
    while (i <= n0) {
      r0 <- rnorm(nSNP, 0, 1)
      r1 <- R1%*%r0 
      r2 <- ifelse(r1<cutoff0, 1, 0)
      
      r0 <- rnorm(nSNP, 0, 1)
      r1 <- R1%*%r0 
      r3 <- ifelse(r1<cutoff0, 1, 0)
      
      r4 <- r2+r3 ## RV
      
      if(nSNP0==0){
        E <- sum(r4)*log(eta)+rnorm(1, 0, 1)
      }else{
        E <- sum(r4[1:nSNP0])*log(eta)+rnorm(1, 0, 1)
      }
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
      r0 <- rnorm(nSNP, 0, 1)
      r1 <- R1%*%r0 
      r2 <- ifelse(r1<cutoff0, 1, 0)
      
      r0 <- rnorm(nSNP, 0, 1)
      r1 <- R1%*%r0 
      r3 <- ifelse(r1<cutoff0, 1, 0)
      
      r4 <- r2+r3 ## RV
      
      if(nSNP0==0){
        E <- sum(r4)*log(eta)+rnorm(1, 0, 1)
      }else{
        E <- sum(r4[1:nSNP0])*log(eta)+rnorm(1, 0, 1)
      }
      ge <- E*r4
      tmp <- exp(alpha0+sum(c(betag, betagx)*c(r4,ge)))
      pr <- tmp/(1+tmp)
      Y1 <- sample(c(0, 1), 1, prob = c(1-pr, pr))
      if(Y1==1){
        G0[i, ] <- r4
        X[i] <- E
        i <- i+1
      }
    }
    
    re <- list(Y=Y, G=G0, X=cbind(X, X1, X2, X3))
    return(re) 
  }
  
  fun <- function(i){
    Data <- simLDRareSNP.joint(para)
    return(Data)
  }
  
  library(parallel)
  system.time(data <- mclapply(1:M, fun, mc.cores=40))
  save(data, file=paste0(dir0, "/", n0, test, "_", set, "_", nSNP, "_", "LD", "_", eta, ".RData"))
  rm(list = ls())
  
}


