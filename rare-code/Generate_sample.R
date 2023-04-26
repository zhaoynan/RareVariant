args=commandArgs(trailingOnly=TRUE)
M = as.numeric(args[1])
n0 = as.numeric(args[2]) # of controls/cases
nSNP = as.numeric(args[3]) # of non-casual SNPs
rho = as.numeric(args[4]) #AR1(rho) between casual SNPs and non-causal SNPs
eta = as.numeric(args[5]) # the degree of correlation between E and causal SNPs
test = args[6] ## null or alter
set = args[7]  ## null--joint/inter  

######################## Generate data ##########################################
library(Matrix)
#variant <- "CV-RV"
variant <- "RV"
  
load(paste0("~/Rare-variant", "/GData", "/",variant, "/Gdata", rho, "_", 1, ".RData"))
dir0 <- paste0("~/Rare-variant", "/", "Joint_test", "/", variant,"/", test)
#dir0 <- paste0("~/Rare-variant", "/", "Inter_test", "/", variant,"/", test)
setwd(dir0)
source("~/Rare-variant/SimAR1rareSNP.R")
 
############## parameters settup #######
set.seed(1234)
nSNP0 <- 8 
if(test=="null"){
  if(set=="joint"){
    ORg <- ORgx <- rep(1, nSNP0)
  }else{
    ORgx <- rep(1, nSNP0)
    ORg <- runif(nSNP0, 0, 2)
  }
}else{
  if(variant=="RV"){
    # interaction test
    #ORgx <- c(1.5, 1.2, 1.2, 1.1, 0.6, 0.8, 0.6, 0.6) 
    #ORg <- runif(nSNP0, 0, 2)
    # joint test
    #ORgx <- runif(nSNP0, 0.5, 1.3) 
    #ORg <- runif(nSNP0, 0, 2)
    # compare joint test and interaction test
    ORgx <- runif(nSNP0, 0.2, 1.6) #runif(nSNP0, 0.25, 1.5)
    ORg <- rep(1, nSNP0)
  }else{
    # interaction test
    #ORgx <- c(1.2, 1.2, 1.2, 1.2, 0.8, 1.1, 1.1, 0.7) 
    #ORg <- runif(nSNP0, 0, 2)
    # joint test  
    # ORgx <- c(0.8, 1.1, 1.1, 1.1, 1.1, 1.1, 0.8, 0.8)
    # ORg <- runif(nSNP0, 0.5, 1.3)
    # compare joint test and interaction test
    ORgx <- c(0.7, 1.1, 1.1, 1.1, 1.1, 1.1, 0.7, 0.7)
    ORg <- rep(1, nSNP0)
  }
}

ORx <- c(0.25, 0.05, 0.5, 0.1)
#################### sample size and prevalence and other parameters settup #########
n1 <- n0
n <- n0+n1
f0 <- 0.01
###################### generate Y and X###########################
YX <- simAR1RareYX(ORg, ORx, eta, ORgx, Gdata, f0)
fn <- function(i){
  id0 <- sample(which(YX$Y==0), n0, replace = FALSE)
  id1 <- sample(which(YX$Y==1), n1, replace = FALSE)
  id <- c(id0, id1)
  
  re <- list(Y=YX$Y[id], G=Gdata[id, ], X=YX$X[id, ])
  return(re)
}

library(parallel)
system.time(data <- mclapply(1:M, fn, mc.cores=30))

if(nSNP > 0){
  if(variant=="CV-RV"){
    MAFslow=0.001; MAFsup=0.45
  }else{
    MAFslow=0.001; MAFsup=0.05
  }
  MAFs <- runif(nSNP, MAFslow, MAFsup)
  ######################## generate noncausal variants of n samples ###############
  para <- list(n=(n*M), nSNP=nSNP, rho=rho, MAFs=MAFs)
  Gnoise <- simAR1RareSNPnoise(para)
  for(i in 1:M){
    data[[i]]$G <- cbind(data[[i]]$G, Gnoise[((i-1)*n+1):(i*n), ])
  }
}
save(data, file=paste0(dir0, "/", n0, test, "_", set, "_", nSNP, "_", rho, "_", eta, ".RData"))
rm(list = ls())

