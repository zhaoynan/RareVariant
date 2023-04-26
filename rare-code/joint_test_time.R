args=commandArgs(trailingOnly=TRUE)
M = as.numeric(args[1]) # number of replications
n0 = as.numeric(args[2])
nSNP = as.numeric(args[3]) # number of non-causal variants 
rho = as.numeric(args[4]) # the parameter for LD
eta = as.numeric(args[5])
test = args[6]
set = args[7]
variant = args[8]


library(gee)
library(MASS, Matrix)
library(mvtnorm, fields)
library(aSPU, magrittr)
library(pryr)
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

######### aGE, IndaGE, and waGE ####################
fn.aGE <- function(x){
  p.aGE <- aGE.joint(data[[x]], pow = c(1:6, Inf), n.perm=2000, 
                     mainweights=weights, interweights=weights)
  return(p.aGE)
} 

mem_aGE <- rep(NA, M)

for(i in 1:M){
  mem_aGE[i] <- mem_change(fn.aGE(i))
}
write.table(mem_aGE, file='mem_aGE.txt', quote=FALSE,append=FALSE)  
rm(list = ls())
rstudioapi::restartSession()
#######################################
fn.IndaGE <- function(x){
  p.IndaGE <- IndaGE.joint(data[[x]], pow = c(1:6, Inf), n.perm=2000, 
                           mainweights=weights, interweights=weights)
  return(p.IndaGE)
}
 
mem_IndaGE <- rep(NA, M)
for(i in 1:M){
  mem_IndaGE[i] <- mem_change(fn.IndaGE(i))
}
write.table(mem_IndaGE, file='mem_IndaGE.txt', quote=FALSE,append=FALSE)
rm(list = ls())
rstudioapi::restartSession()
#########################################
fn.aGE <- function(x){
  p.aGE <- aGE.joint(data[[x]], pow = c(1:6, Inf), n.perm=2000, 
                     mainweights=weights, interweights=weights)
  return(p.aGE)
} 

fn.IndaGE <- function(x){
  p.IndaGE <- IndaGE.joint(data[[x]], pow = c(1:6, Inf), n.perm=2000, 
                           mainweights=weights, interweights=weights)
  return(p.IndaGE)
}
p_aGE <- p_IndaGE <- matrix(NA, ncol=10, nrow = M)
for(i in 1:M){
  p_aGE[i, ] <- fn.aGE(i)
  p_IndaGE[i, ] <- fn.IndaGE(i)
}

write.table(p_aGE, file='p_aGE.txt', quote=FALSE,append=FALSE)
write.table(p_IndaGE, file='p_IndaGE.txt', quote=FALSE,append=FALSE)
rm(list = ls())
rstudioapi::restartSession()
##############################################
fn.waGE <- function(x){
  p.waGE <- WaGE.joint(data[[x]], p_aGE[x,], p_IndaGE[x,], n.perm=2000, alpha=0.1)
  return(p.waGE)
}

mem_waGE <- rep(NA, M)

p_aGE <- read.csv('p_aGE.txt' , header = TRUE, sep = "")
p_IndaGE <- read.csv('p_IndaGE.txt', header = TRUE, sep = "")

for(i in 1:M){
  mem_waGE[i] <-  mem_change(fn.waGE(i))
}

write.table(mem_waGE, file='mem_waGE.txt', quote=FALSE,append=FALSE)
rm(list = ls())
rstudioapi::restartSession()
########### rareGE, IndrareGE, and wrareGE ##############

fn.rareGE <- function(x){
  p.rareGE <- rareGE.joint(data[[x]], n.perm=2000, 
                           mainweights=weights, interweights=weights)
  return(p.rareGE)
}
mem_rareGE <- rep(NA, M) 

for(i in 1:M){
  mem_rareGE[i] <-  mem_change(fn.rareGE(i))
}

write.table(mem_rareGE, file='mem_rareGE.txt', quote=FALSE,append=FALSE)
rm(list = ls())
rstudioapi::restartSession()
#####################################
fn.IndrareGE <- function(x){
  p.IndrareGE <- IndrareGE.joint(data[[x]], n.perm =2000, 
                                 mainweights=weights, interweights=weights)
  return(p.IndrareGE)
}

mem_IndrareGE <- rep(NA, M) 
for(i in 1:M){
  mem_IndrareGE[i] <-  mem_change(fn.IndrareGE(i))
}
write.table(mem_IndrareGE, file='mem_IndrareGE.txt', quote=FALSE,append=FALSE)

rm(list = ls())
rstudioapi::restartSession()
######################################
fn.rareGE <- function(x){
  p.rareGE <- rareGE.joint(data[[x]], n.perm=2000, 
                           mainweights=weights, interweights=weights)
  return(p.rareGE)
}

fn.IndrareGE <- function(x){
  p.IndrareGE <- IndrareGE.joint(data[[x]], n.perm =2000, 
                                 mainweights=weights, interweights=weights)
  return(p.IndrareGE)
}
p_rareGE <- p_IndrareGE <- rep(NA, M) 
for(i in 1:M){
  p_rareGE[i] <- fn.rareGE(i)
  p_IndrareGE[i] <- fn.IndrareGE(i)
}
write.table(p_rareGE, file='p_rareGE.txt', quote=FALSE,append=FALSE)
write.table(p_IndrareGE, file='p_IndrareGE.txt', quote=FALSE,append=FALSE)
rm(list = ls())
rstudioapi::restartSession()
########################################
p_rareGE <- as.vector(read.csv('p_rareGE.txt' , header = TRUE, sep = ""))$x 
p_IndrareGE <- as.vector(read.csv('p_IndrareGE.txt', header = TRUE, sep = ""))$x

fn.wrareGE <- function(x){
  p.wrareGE <- WrareGE.joint(data[[x]], p_rareGE[x], p_IndrareGE[x], n.perm=2000, alpha=0.1)
  return(p.wrareGE)
}
mem_wrareGE <- rep(NA, M)

for(i in 1:M){
  mem_wrareGE[i] <- mem_change(fn.wrareGE(i))
}
write.table(mem_wrareGE, file='mem_wrareGE.txt', quote=FALSE,append=FALSE)

rm(list = ls())
rstudioapi::restartSession()

################################################################
p_aGE <- p_IndaGE <- p_waGE <- matrix(NA, ncol=10, nrow = M)
p_rareGE <- p_IndrareGE <- p_wrareGE <- rep(NA, M)
#time <- Ind.time <- W.time <- matrix(NA, ncol = 2, nrow = M)
mem <- Ind.mem <- W.mem <- matrix(NA, ncol = 2, nrow = M)


for(i in 1:M){
  sta.aGE <- Sys.time()
  p_aGE[i, ] <- fn.aGE(i)
  end.aGE <- Sys.time()
  
  sta.rareGE <- Sys.time()
  p_rareGE[i] <- fn.rareGE(i)
  end.rareGE <- Sys.time()
  #############################
  time[i,  ] <- c(end.aGE-sta.aGE, end.rareGE-sta.rareGE)
  ############################# 
  sta.IndaGE <- Sys.time()
  p_IndaGE[i, ] <- fn.IndaGE(i)
  end.IndaGE <- Sys.time()
  
  sta.IndrareGE <- Sys.time()
  p_IndrareGE[i] <- fn.IndrareGE(i)
  end.IndrareGE <- Sys.time()
  #############################
  Ind.time[i,  ] <- c(end.IndaGE-sta.IndaGE, end.IndrareGE-sta.IndrareGE) 
  ############################# 
  sta.waGE <- Sys.time()
  p_waGE[i, ] <- fn.waGE(i)
  end.waGE <- Sys.time()
  
  sta.wrareGE <- Sys.time()
  p_wrareGE[i] <- fn.wrareGE(i)
  end.wrareGE <- Sys.time()
  #############################
  W.time[i,  ] <- c(end.waGE-sta.waGE, end.wrareGE-sta.wrareGE)  
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
 
name <- c(c(paste("aGEj", pow = c(1:6, Inf), sep = ""), "aGEj", "aGEj_minp", "aGEj_fisher"), 
          c(paste("IndaGEj", pow = c(1:6, Inf), sep = ""), "IndaGEj","IndaGEj_minp", "IndaGEj_fisher"), 
          c(paste("waGEj", pow = c(1:6, Inf), sep = ""), "waGEj","waGEj_minp", "waGEj_fisher"), 
          c("rareGE", "IndrareGE", "wrareGE"))
colnames(p_value) <- name
rownames(power) <- name

time.inter <- cbind(time, Ind.time, W.time)
mean.time <- colMeans( time.inter, na.rm = T)
se.time <- apply(time.inter, 2, sd, na.rm=TRUE)
re.time <- rbind(mean.time, se.time) 
colnames(re.time) <- c("aGE", "rareGE", "IndaGE", "IndrareGE", "waGE", "wrareGE") 

write.table(power, file='power_joint_test.txt',  quote=FALSE,append=FALSE)
write.table(p_value, file='Pvalue_joint_test.txt', quote=FALSE,append=FALSE) ### For ROC and AUC
write.table(time.inter, file='time_joint.txt', quote=FALSE,append=FALSE) 
write.table(re.time, file='re_time.txt', quote=FALSE,append=FALSE) 

rm(list = ls())
