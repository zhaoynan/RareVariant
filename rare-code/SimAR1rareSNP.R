########### simulate a haplotype 
simHaplotype <- function(nSNP, R, MAFs){
  #nSNP: the number of SNPs in the simulated haplotype;
  #R: the relationship matrix of SNPs in the simulated haplotype;
  #MAFs: the MAF of the nSNP SNPs;
  if(ncol(R) != nrow(R)) stop("R should be a square matrix!")
  if(nSNP != ncol(R)) stop("The dimension of R is the same as the number of SNPs!")
  
  eigen.R <- eigen(R, symmetric = T)
  R1 <- eigen.R$vectors%*% diag(sqrt(eigen.R$values))
  X0 <- rnorm(nSNP, 0, 1) #: X0 ~ MVN(0, I)
  X1 <- R1 %*% X0   #: X1 ~ MVN(0, R)
  
  cutoff <- qnorm(MAFs)
  X2 <- ifelse(X1 < cutoff, 1, 0)  
  return(X2)
}
 
#### Simulate causal SNPs for an individual:SNPs from a latent multivariate Gaussian variable with
####an AR1(rho) correlation structure;
simAR1RareSNP0 <- function(para){
  
  # nSNP0: # of causal SNPs; nSNP: # of non-causal SNPs
  # rho0: the parameter in AR1(rho) for the latent multivariate Gaussian variable to be discretized to causal SNPs;
  #       rho0=0 means all SNPs are independent;
  # rho1: the parameter in AR1(rho) for the latent multivariate Gaussian variable to be discretized to non-causal SNPs;
  #       rho=0 means all SNPs are independent; 
  # MAF0slow, MAF0sup: MAFs for the causal SNPs from  Unif(MAF0slow, MAF0sup);
  # MAFslow, MAFsup: MAFs of the non-causal SNPs are drawn from Unif(MAFslow, MAFsup);
  
  nSNP0 <- para$nSNP0; 
  rho0 <- para$rho0;  
  MAF0s <- para$MAF0s; 
 
  N <- para$N         # the number of pupulation
  
  R0 <- matrix(1, nrow=nSNP0, ncol=nSNP0)
  for(i in 1:nSNP0)
    for(j in 1:nSNP0)
      if (i!=j) R0[i, j]<- rho0^(abs(i-j))
  
  fun0 <- function(i){
    data <- simHaplotype(nSNP0, R0, MAF0s)
    return(data)
  } 
  
  library(parallel)
  system.time(Hap0 <- mclapply(1:N, fun0, mc.cores=25))
  G00 <- Matrix(unlist(Hap0), nrow = N, byrow = T)
  rm(Hap0)
  system.time(Hap1 <- mclapply(1:N, fun0, mc.cores=25))
  G01 <- Matrix(unlist(Hap1), nrow = N, byrow = T)
  G0 <- G00+G01
  rm(Hap1) 
  return(G0)
}
############################################
# generate non-causal causal SNPs, the block of causal SNPs and the block of non-causal ones are independent. 
simAR1RareSNPnoise <- function(para){
  
  # nSNP0: # of causal SNPs; nSNP: # of non-causal SNPs
  # rho0: the parameter in AR1(rho) for the latent multivariate Gaussian variable to be discretized to causal SNPs;
  #       rho0=0 means all SNPs are independent;
  # rho1: the parameter in AR1(rho) for the latent multivariate Gaussian variable to be discretized to non-causal SNPs;
  #       rho=0 means all SNPs are independent; 
  # MAF0slow, MAF0sup: MAFs for the causal SNPs from  Unif(MAF0slow, MAF0sup);
  # MAFslow, MAFsup: MAFs of the non-causal SNPs are drawn from Unif(MAFslow, MAFsup);
  
  nSNP <- para$nSNP; 
  rho <- para$rho;  
  MAFs <- para$MAFs; 
  n <- para$n         # the number of pupulation
  
  R <- matrix(1, nrow=nSNP, ncol=nSNP)
  for(i in 1:nSNP)
    for(j in 1:nSNP)
      if (i!=j) R[i, j]<- rho^(abs(i-j))
  
  fun <- function(i){
    data <- simHaplotype(nSNP, R, MAFs)
    return(data)
  } 
  
  library(parallel)
  system.time(Hap0 <- mclapply(1:n, fun, mc.cores=30))
  G00 <- Matrix(unlist(Hap0), nrow = n, byrow = T)
  rm(Hap0)
  system.time(Hap1 <- mclapply(1:n, fun, mc.cores=30))
  G01 <- Matrix(unlist(Hap1), nrow = n, byrow = T)
  G0 <- G00+G01
  rm(Hap1) 
  return(G0)
}

### generate phenotype, covariate ###################
simAR1RareYX <- function(ORg, ORx, eta, ORgx, G0, f0){
  # ORg:association in OR between the causal SNPs and outcome; implicitly,  # of causal SNPs = length(ORg)
  # ORx:association in OR between the covariates and outcome; implicitly,  # of covariates = length(ORx);
  # the first value is environmental factor
  # ORgx: association in OR between the G-E interactions and outcome; implicitly
  # eta: association between E and causal SNPs
  # G: genotype data for N individuals, and the first nSNP0 are causal SNPs
  #nGE: the number of SNPs related to E
  nSNP0 <- length(ORg)
  nCov <- length(ORx)
  N <- nrow(G0)
  
  betag <- log(ORg)
  betax <- ORx
  betagx <- log(ORgx)
  beta <- c(betax, betag, betagx)
 
  X1 <- rnorm(N, mean = 30, sd=5.7)
  X2 <- rbinom(N, 1, prob = 0.5)
  X3 <- rnorm(N, mean = 160, sd=9.4)
  
  E <- as.vector(scale(G0, center=T, scale=F)%*%rep(log(eta), nSNP0)+rnorm(N, 0, 1))
  X4 <- cbind(E, X1, X2, X3)
  Z <- cbind(X4, G0, E*G0)
  
  fn <- function(alpha0){
    tmp <- exp(alpha0 +as.vector(Z%*%beta))
    p <- tmp/(1+tmp)
    return(mean(p)-f0)
  }
  alpha0 <- uniroot(fn, c(-100, 100))$root
  tmp <- exp(alpha0 + as.vector(Z%*%beta))
  mu <- tmp/(1+tmp)
  Y <- rbinom(N, 1, prob=mu) 
  
  re <- list(Y=Y, X=X4, alpha0=alpha0)
  return(re)
}

####
simLDRareSNP <- function(para){
  nSNP0 <- para$nSNP; 
  R0 <- para$R;  
  MAF0s <- para$MAF0s; 
  N <- para$N         # the number of pupulation
 
  fun0 <- function(i){
    data <- simHaplotype(nSNP0, R0, MAF0s)
    return(data)
  } 
  
  library(parallel)
  system.time(Hap0 <- mclapply(1:N, fun0, mc.cores=30))
  G00 <- Matrix(unlist(Hap0), nrow = N, byrow = T)
  rm(Hap0)
  system.time(Hap1 <- mclapply(1:N, fun0, mc.cores=30))
  G01 <- Matrix(unlist(Hap1), nrow = N, byrow = T)
  G0 <- G00+G01
  rm(Hap1) 
  return(G0)
}

