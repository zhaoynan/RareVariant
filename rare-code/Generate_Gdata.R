# Generate the csRVs of N objectives
 
rho <- 0.9
eta <- 1
 
library(Matrix)
source("~/Rare-variant/SimAR1rareSNP.R")
set.seed(1234)
N <- 1e6
nSNP0 <- 8
MAF0slow=0.001; MAF0sup=0.45
 
MAF0s <- runif(nSNP0, MAF0slow, MAF0sup)
para <- list(N=N, nSNP0=nSNP0, rho0=rho,  MAF0s=MAF0s)
 
G <- list()
for(i in 1:10){
  system.time(G[[i]] <- simAR1RareSNP0(para))
}
Gdata <- do.call("rbind", G)
rm(G)
save(Gdata, file=paste0("~/Rare-variant", "/GData",  "/CV-RV","/Gdata", rho, "_", eta, ".RData"))
#save(Gdata, file=paste0("~/Rare-variant", "/GData",  "/RV","/Gdata", rho, "_", eta, ".RData"))
rm(list = ls())

