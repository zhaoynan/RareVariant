
library(data.table)
Pheno_covariate_data <- fread(file = "~/PretermData/PLINK_filtered_dataset/PD_QC_MCpairs_phenotype.txt", header = TRUE, sep = " ")  ## 1626 samples have maternal information
BMINA_ID <- which(is.na(Pheno_covariate_data$maternal_bmi))
pheno_deleteNA <- Pheno_covariate_data[-BMINA_ID, ]

phenotype <- pheno_deleteNA$premie; table(phenotype) ## case: 665, control: 855
BMI <- pheno_deleteNA$maternal_bmi; m_height <- pheno_deleteNA$maternal_height;  m_age <- pheno_deleteNA$maternal_age


#Raw_data <- fread(file = "~/PretermData/PLINK_filtered_dataset/PD_QC_MCpairs_additive.raw", header = TRUE, sep = " ")   # The SNP information of all individual, including  offspring and maternal genotypes
#n_m <- n_c <- nrow(Raw_data)/2
#MGdata <- Raw_data[1:n_m, ]
#rm(Raw_data)
#m_SNP <- gsub('.{2}$', '', colnames(MGdata)[-(1:6)])
#colnames(MGdata)[-(1:6)] <- m_SNP
#MGdata <- as.matrix(MGdata[, -(1:6)])
#save(MGdata, file="~/PretermData/PLINK_filtered_dataset/PD_MGdata.RData")

load(file="~/PretermData/PLINK_filtered_dataset/PD_MGdata.RData")

PD_QC_MCpairs_pheno_BMI.bim <- read.table(file="~/PretermData/PLINK_filtered_dataset/PD_QC_MCpairs_pheno_BMI.bim")
all_SNP <- PD_QC_MCpairs_pheno_BMI.bim$V2
chro <- PD_QC_MCpairs_pheno_BMI.bim$V1
position <- PD_QC_MCpairs_pheno_BMI.bim$V4 

rm(Pheno_covariate_data, PD_QC_MCpairs_pheno_BMI.bim, pheno_deleteNA, all_SNP)
rm(position, BMINA_ID)

# check if covariates NA
all(is.na(m_height) == FALSE) ## true
all(is.na(m_age) == FALSE) ## false so complete NA using mean values
m_age[is.na(m_age)] <- mean(m_age, na.rm = TRUE)

source("~/Rare-variant/Mixed_chi_squared_pvalue.R")
source("~/Rare-variant/Statistic_functions.R")
#source("~/Rare-variant/aGE_functions.R")

dir <- "~/Rare-variant/Real_data_result"

# calculate the number of windows
L <- rep(NA, 22); dis <- 30; m <- 60
for(j in 1:22){
  SNPs <- which(chro == j)
  L[j] <- floor((length(SNPs)-m)/dis+1) 
}
num_M <- sum(L) ## 18128; 0.05/num_M=2.758164e-06

####  joint test ############
for(j in 1:22){
  ## Whole chromosome analysis
  SNPs <- which(chro == j)
  
  re <- list()
  #dis <- 10; m <- 30
  fn <- function(x){
    id_SNP <- SNPs[(1+(x-1)*dis):(m+(x-1)*dis)]
    re$X <- cbind(BMI, m_height, m_age)
    re$G <- as.matrix(MGdata)[, id_SNP]
    re$Y <- phenotype
    
    pval <- Get_joint_realdata(re, alpha=0.1,pow = c(1:6, Inf), n.perm =1000,
                               mainweights=poweights, interweights=poweights)
    p.aGE <- pval$aGE[1:10]
    p.raGE <- pval$aGE[11:20]
    p.aGE.EB <- pval$aGE[21:30]
    p.rareGE <- pval$rareGE
    rho <- pval$rho
    p0 <- pval$p0
    re <- list(p.aGE=p.aGE, p.raGE=p.raGE, p.aGE.EB=p.aGE.EB ,p.rareGE=p.rareGE, rho=rho, p0=p0)
    return(re)
  }
  
  M <- floor((length(SNPs)-m)/dis+1) 
  
  finally.result <- parallel::mclapply(1:M, fn, mc.cores=getOption("mc.cores", 25L))
  
  p_aGE <- p_raGE <- p_aGE.EB <- matrix(NA, ncol = 10, nrow=M)
  p_rareGE <- matrix(NA, ncol = 3, nrow=M)
  rho <- pind <- rep(NA, M)
  for(i in 1:M){
    p_aGE[i, ] <- finally.result[[i]]$p.aGE
    p_raGE[i, ] <- finally.result[[i]]$p.raGE
    p_aGE.EB[i, ] <- finally.result[[i]]$p.aGE.EB
    p_rareGE[i, ] <- finally.result[[i]]$p.rareGE
    rho[i] <- finally.result[[i]]$rho
    pind[i] <- finally.result[[i]]$p0
  }
  
  id_SNP <- SNPs[(length(SNPs)-m):length(SNPs)]
  re$X <- cbind(BMI, m_height, m_age)
  re$G <- as.matrix(MGdata)[, id_SNP]
  re$Y <- phenotype
  p_finally <- Get_joint_realdata(re, alpha=0.1,pow = c(1:6, Inf), n.perm =1000, 
                                  mainweights=poweights, interweights=poweights) 
  p.aGE <- p_finally$aGE[1:10]
  p.raGE <- p_finally$aGE[11:20]
  p.aGE.EB <- p_finally$aGE[21:30]
  p.rareGE <- p_finally$rareGE
  
  p_value <- cbind(rbind(p_aGE, p.aGE), rbind(p_raGE, p.raGE), rbind(p_aGE.EB, p.aGE.EB), rbind(p_rareGE, p.rareGE))
  rho <- c(rho, p_finally$rho)
  pind <- c(pind, p_finally$p0)
  
  write.table(p_value, file = paste0(dir, "/joint_result", j, ".txt"), quote = FALSE,row.names = FALSE, append=FALSE)
  write.table(cbind(rho, pind), file = paste0(dir, "/joint_rho", j, ".txt"), quote = FALSE, row.names = FALSE, append=FALSE)
}

### inter test ############
for(j in 1:22){
  ## Whole chromosome analysis
  SNPs <- which(chro == j)
  
  re <- list()
  #dis <- 10; m <- 30
  finally <- function(x){
    id_SNP <- SNPs[(1+(x-1)*dis):(m+(x-1)*dis)]
    re$X <- cbind(BMI, m_height, m_age)
    re$G <- as.matrix(MGdata)[, id_SNP]
    re$Y <- phenotype
    
    pval <- Get_inter_EB_glmm(re, alpha=0.1,pow = c(1:6, Inf), n.perm =1000, 
                              mainweights=poweights, interweights=poweights)
    p.aGE <- pval$aGE[1:9]
    p.raGE <- pval$aGE[10:18]
    p.aGE.EB <- pval$aGE[19:27]
    p.rareGE <- pval$rareGE
    rho <- pval$rho
    p0 <- pval$p0
    phi <- pval$phi
    re <- list(p.aGE=p.aGE, p.raGE=p.raGE, p.aGE.EB=p.aGE.EB ,p.rareGE=p.rareGE,
               rho=rho, p0=p0, phi=phi)
    return(re)
  }
  
  M <- floor((length(SNPs)-m)/dis+1) 
  
  finally.result <- parallel::mclapply(1:M, finally, mc.cores=getOption("mc.cores", 25L))
  
  p_aGE <- p_raGE <- p_aGE.EB <- matrix(NA, ncol = 9, nrow=M)
  p_rareGE <- matrix(NA, ncol = 3, nrow=M)
  rho <- pind <- Phi <- rep(NA, M)
  for(i in 1:M){
    p_aGE[i, ] <- finally.result[[i]]$p.aGE
    p_raGE[i, ] <- finally.result[[i]]$p.raGE
    p_aGE.EB[i, ] <- finally.result[[i]]$p.aGE.EB
    p_rareGE[i, ] <- finally.result[[i]]$p.rareGE[1:3]
    rho[i] <- finally.result[[i]]$rho
    pind[i] <- finally.result[[i]]$p0
    Phi[i] <- finally.result[[i]]$phi
  }
  
  id_SNP <- SNPs[(length(SNPs)-m):length(SNPs)]
  re$X <- cbind(BMI, m_height, m_age)
  re$G <- as.matrix(MGdata)[, id_SNP]
  re$Y <- phenotype
  p_finally <- Get_inter_EB_glmm(re, alpha=0.1,pow = c(1:6, Inf), n.perm =1000,
                                  mainweights=poweights, interweights=poweights) 
  p.aGE <- p_finally$aGE[1:9]
  p.raGE <- p_finally$aGE[10:18]
  p.aGE.EB <- p_finally$aGE[19:27]
  p.rareGE <- p_finally$rareGE[1:3]
  
  p_value <- cbind(rbind(p_aGE, p.aGE), rbind(p_raGE, p.raGE), rbind(p_aGE.EB, p.aGE.EB), rbind(p_rareGE, p.rareGE))
  rho <- c(rho, p_finally$rho)
  pind <- c(pind, p_finally$p0)
  Phi <- c(Phi, p_finally$phi)
  write.table(p_value, file = paste0(dir, "/inter_result", j, ".txt"), quote = FALSE,row.names = FALSE, append=FALSE)
  write.table(cbind(rho, pind, Phi), file = paste0(dir, "/inter_rho", j, ".txt"), quote = FALSE, row.names = FALSE, append=FALSE)
}

######### main test ###########
for(j in 1:22){
  ## Whole chromosome analysis
  SNPs <- which(chro == j)
  
  re <- list()
  #dis <- 10; m <- 30
  finally <- function(x){
    id_SNP <- SNPs[(1+(x-1)*dis):(m+(x-1)*dis)]
    re$X <- cbind(BMI, m_height, m_age)
    re$G <- as.matrix(MGdata)[, id_SNP]
    re$Y <- phenotype
    
    pval <- Get_main_realdata(re, alpha=0.1,pow = c(1:8, Inf), n.perm =2000, 
                               mainweights=poweights, interweights=poweights)
    p.aspu <- pval$aspu[1:11]
    p.raspu <- pval$aspu[12:22]
    p.aspu.EB <- pval$aspu[23:33]
    rho <- pval$rho
    p0 <- pval$p0
    re <- list(p.aspu=p.aspu, p.raspu=p.raspu, p.aspu.EB=p.aspu.EB, rho=rho, p0=p0)
    return(re)
  }
  
  M <- floor((length(SNPs)-m)/dis+1) 
  
  finally.result <- parallel::mclapply(1:M, finally, mc.cores=getOption("mc.cores", 20L))
  
  p_aspu <- p_raspu <- p_aspu.EB <- matrix(NA, ncol = 11, nrow=M)
  rho <- pind <- rep(NA, M)
  for(i in 1:M){
    p_aspu[i, ] <- finally.result[[i]]$p.aspu
    p_raspu[i, ] <- finally.result[[i]]$p.raspu
    p_aspu.EB[i, ] <- finally.result[[i]]$p.aspu.EB

    rho[i] <- finally.result[[i]]$rho
    pind[i] <- finally.result[[i]]$p0
  }
  
  id_SNP <- SNPs[(length(SNPs)-m):length(SNPs)]
  re$X <- cbind(BMI, m_height, m_age)
  re$G <- as.matrix(MGdata)[, id_SNP]
  re$Y <- phenotype
  p_finally <- Get_main_realdata(re, alpha=0.1,pow = c(1:8, Inf), n.perm =2000, 
                                  mainweights=poweights, interweights=poweights) 
  p.aspu <- p_finally$aspu[1:11]
  p.raspu <- p_finally$aspu[12:22]
  p.aspu.EB <- p_finally$aspu[23:33]
  
  p_value <- cbind(rbind(p_aspu, p.aspu), rbind(p_raspu, p.raspu), rbind(p_aspu.EB, p.aspu.EB))
  rho <- c(rho, p_finally$rho)
  pind <- c(pind, p_finally$p0)
  
  write.table(p_value, file = paste0(dir, "/Main","/main_result", j, ".txt"), quote = FALSE,row.names = FALSE, append=FALSE)
  write.table(cbind(rho, pind), file = paste0(dir, "/Main","/main_rho", j, ".txt"), quote = FALSE, row.names = FALSE, append=FALSE)
}
