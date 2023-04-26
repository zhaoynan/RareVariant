args=commandArgs(trailingOnly=TRUE)
M = as.numeric(args[1]) # number of replications
nSNP = as.numeric(args[2]) # number of non-causal variants 
rho = as.numeric(args[3]) # the parameter for LD
eta = as.numeric(args[4])
test = args[5]
set = args[6]

library(aSPU)
dir0 <- paste0("~/Rare-variant", "/", "Joint_test")
setwd(dir0)
set.seed(1234)

source("~/Rare-variant/Mixed_chi_squared_pvalue.R")
source("~/Rare-variant/HC_pvalue.R")
source("~/Rare-variant/Statistic_functions.R")
source("~/Rare-variant/aGE_functions.R")

dir <- paste0(dir0, "/", n0, test, "_", set, "_", nSNP, "_", rho, "_", eta, ".RData")

dir <- paste0(dir0, "/", set,  nSNP, "_", rho, "_", eta, "_", test)
dir.create(dir, showWarnings = FALSE)
setwd(dir)
load(paste0(dir0,"/", test, "/", n0, test, "_", set, "_", nSNP, "_", rho, "_", eta, ".RData"))

finally <- function(x){
  n <- dim(data[[x]]$G)[1]
  m <- dim(data[[x]]$G)[2]
  k <- 10
  stat_score.glm <- Get_joint_HC_glm(data[[x]] ,alpha=0.05)
  rho <- stat_score.glm$rho
  rho.roc <- stat_score.glm$rho.roc
  
  Sigma.Z.s <- stat_score.glm$Sigma$Sigma.Z.s
  Sigma.Z.star <- stat_score.glm$Sigma$Sigma.Z.s.star
  
  p.s <- stat_score.glm$p_score$p.s 
  p.star <- stat_score.glm$p_score$p.s.star
  
  stat.s <- stat_WRHC(Sigma.Z.s, p.s, k)
  stat.star <- stat_WRHC(Sigma.Z.star, p.star, k)
  
  t.s <- Get_tk(data[[x]], Sigma.Z.s, stat.s, k)
  t.star <- Get_tk(data[[x]], Sigma.Z.star, stat.star, k)
  
  p.s <- Get_pvalue(Sigma.Z.s, t.s, k)$p.value
  p.star <- Get_pvalue(Sigma.Z.star, t.star, k)$p.value
  p.EB <- ifelse(rho==0, p.star, p.s)
  p.roc <- ifelse(rho.roc==0, p.star, p.s)
  
  p.HC <- list(p.s=p.s, p.star=p.star, p.EB=p.EB, p.roc=p.roc)
  return(p.HC)
}

library(parallel)
system.time(finally.result <- mclapply(1:M, finally, mc.cores=10))

p_value <- matrix(NA, ncol = 16, nrow=M)
rho <- rho.roc <- numeric(M)

for(i in 1:M){
  p_value[i, ] <- c(finally.result[[i]]$p.st, finally.result[[i]]$p.bt, finally.result[[i]]$p.fisher, finally.result[[i]]$p.min)
  rho[i] <- finally.result[[i]]$rho
  rho.roc[i] <- finally.result[[i]]$rho.roc
}
power_01 <- colMeans(p_value <= 0.01, na.rm=TRUE)
power_03 <- colMeans(p_value <= 0.03, na.rm=TRUE)
power_05 <- colMeans(p_value <= 0.05, na.rm=TRUE)
power_07 <- colMeans(p_value <= 0.07, na.rm=TRUE)
power_1 <- colMeans(p_value <= 0.1, na.rm=TRUE)
power_3 <- colMeans(p_value <= 0.3, na.rm=TRUE)
power_5 <- colMeans(p_value <= 0.5, na.rm=TRUE)
power_7 <- colMeans(p_value <= 0.7, na.rm=TRUE)
power_9 <- colMeans(p_value <= 0.9, na.rm=TRUE)
power.ind <- mean(rho)
power.ROC <- mean(rho.roc)
power <- cbind(power_01, power_03,power_05, power_07,power_1, power_3, power_5, power_7, power_9)
write.table(power, file='Joint.txt',row.names=c("ST", "RST", "EBST","ROCST","BT", "RBT", "EBBT",
                                                "ROCBT","aGEfisher", "RaGEfisher", "EBaGEfisher", "ROCaGEfisher",
                                                "aGEmin", "RaGEmin", "EBaGEmin", "ROCaGEmin"), quote=FALSE,append=FALSE)
write.table(c(power.ind, power.ROC), file='Ind_alltest.txt', quote=FALSE,append=FALSE)
rm(list = ls())


