args=commandArgs(trailingOnly=TRUE)
M = as.numeric(args[1]) # number of replications
n0 = as.numeric(args[2])
nSNP = as.numeric(args[3]) # number of non-causal variants 
rho = as.numeric(args[4]) # the parameter for LD
eta = as.numeric(args[5])
test = args[6]
set = args[7]
variant <- args[8]

library(gee)
library(MASS, Matrix)
library(mvtnorm, fields)
library(aSPU, magrittr)

dir0 <- paste0("~/Rare-variant", "/", "Joint_test","/", variant)
setwd(dir0)
set.seed(1234)

#source("~/Rare-variant/Mixed_chi_squared_pvalue.R")
source("~/Rare-variant/aGE_functions.R")
source("~/Rare-variant/Statistic_functions1.R")
source("~/Rare-variant/Statistic_functions2.R")

 
dir <- paste0(dir0, "/", set,  nSNP, "_", rho, "_", eta, "_", test, n0)
dir.create(dir, showWarnings = FALSE)
setwd(dir)
load(paste0(dir0,"/", test, "/", n0, test, "_", set, "_", nSNP, "_", rho, "_", eta, ".RData"))

weights <- ifelse(variant=="RV", wuweights, poweights)

fn.result <- function(x){
  p <- Get_joint_glm(data[[x]], alpha=0.1, pow = c(1:6, Inf), n.perm =2000, 
                      mainweights=weights, interweights=weights) 
    
  return(p)
}

pval.result <- parallel::mclapply(1:M, fn.result, mc.cores=getOption("mc.cores", 40L))

################################################################
p_aGE <- p_IndaGE <- p_waGE <- matrix(NA, ncol=10, nrow = M)
p_rareGE <- p_IndrareGE <- p_wrareGE <- rep(NA, M)
for(i in 1:M){
  p_aGE[i, ] <- pval.result[[i]]$p.aGE
  p_IndaGE[i, ] <- pval.result[[i]]$p.IndaGE
  p_waGE[i, ] <- pval.result[[i]]$p.waGE
  
  p_rareGE[i] <- pval.result[[i]]$p.rareGE
  p_IndrareGE[i] <- pval.result[[i]]$p.IndrareGE
  p_wrareGE[i] <- pval.result[[i]]$p.wrareGE
}
p_value <- cbind(p_aGE, p_IndaGE, p_waGE, p_rareGE, p_IndrareGE, p_wrareGE)
colnames(p_value) <- NULL
power_01 <- colMeans(p_value <= 0.01, na.rm=TRUE)
power_03 <- colMeans(p_value <= 0.03, na.rm=TRUE)
power_05 <- colMeans(p_value <= 0.05, na.rm=TRUE)
power_07 <- colMeans(p_value <= 0.07, na.rm=TRUE)
power_1 <- colMeans(p_value <= 0.1, na.rm=TRUE)
power_3 <- colMeans(p_value <= 0.3, na.rm=TRUE)
power_5 <- colMeans(p_value <= 0.5, na.rm=TRUE)
power_7 <- colMeans(p_value <= 0.7, na.rm=TRUE)
power_9 <- colMeans(p_value <= 0.9, na.rm=TRUE)

power <- cbind(power_01, power_03,power_05, power_07,power_1, power_3, power_5, power_7, power_9)
#time <- cbind(time.aGE.joint, time.IndaGE.joint, time.waGE.joint, time.rareGE.joint, time.IndrareGE.joint, time.wrareGE.joint)

name <- c(c(paste("aGEj", pow = c(1:6, Inf), sep = ""), "aGEj", "aGEj_minp", "aGEj_fisher"), 
          c(paste("IndaGEj", pow = c(1:6, Inf), sep = ""), "IndaGEj","IndaGEj_minp", "IndaGEj_fisher"), 
          c(paste("waGEj", pow = c(1:6, Inf), sep = ""), "waGEj","waGEj_minp", "waGEj_fisher"), 
          c("rareGE", "IndrareGE", "wrareGE"))
colnames(p_value) <- name
rownames(power) <- name
write.table(power, file='power_joint_test.txt',  quote=FALSE,append=FALSE)
write.table(p_value, file='Pvalue_joint_test.txt', quote=FALSE,append=FALSE) ### For ROC and AUC
#write.table(time, file='time_joint_test.txt', quote=FALSE,append=FALSE) 

rm(list = ls())



