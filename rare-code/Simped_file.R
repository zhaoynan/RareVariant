####### To create two files required for Simped to simulate genotype data ###############

######################## To create the pedifile file with N=1e7

## Column 1: The family ID
## Column 2: The individual's ID
#Column 3: the ID of the individual’s father
#Column 4: the ID of the individual’s mother
#Column 5: sex (1=male, 2=female)
dir <- "~/Rare-variant/GData"
n <- 16 # the number of members in a pedigree
FamID <- rep(1, n) 
IID <- seq(1, n, 1)
FatID <- rep(0, n) 
MomID <- rep(0, n) 
sex <- c(2, 1, 1, 2, 2,2, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2)
pedin <- data.frame(FamID=FamID,IID = IID, FatID=FatID, MomID=MomID, sex=sex, G=rep(1, n))
write.table(pedin, file=paste0(dir, "/pedin_unrelated", ".pre"), quote=F, row.names = F, col.names = F)
#ped <- pedigree(IID,FatID, MomID, sex,famid=FamID)
#ped <- ped['1']
#plot(ped)  ## pedigree chart
#related_ped <- with(pedin, pedigree(IID, FatID, MomID, sex))
#Phi <- 2*kinship(related_ped) ## kinship matrix
#########################################################################
######################## To create the parameter file for phenotype data #########
rho <- 0.7
nSNP0 <- 8
names <- paste("pedin_unrelated.pre",  paste0("pedfile_unrelated","_", rho, "_", nSNP0, ".pre")) ### name of pedigree file, name of output file
seed <- paste(23221, 1601, 21001)             ## three random seeds
traid <- 0   # number of columns for affection status/quantitative trait
n.fam <- 625000  # number of replicates/families
n.snp <- paste(nSNP0, 1)    # Total number of marker loci, # of times pattern to be repeated 
re <- 1  # “1” recomb fraction, “2” Kosambi map distance & “3” Haldane map distance  
re.num <- matrix(c(nSNP0-1, rep(0, nSNP0-1)), nrow=1)  # no recombination
## Generate the phenotype
maf <- runif(nSNP0, 0.001, 0.05)
sigma <- (1-rho)*diag(1, nSNP0, nSNP0)+rho*matrix(1, nSNP0, nSNP0)
upper <- sapply(rep(Inf, nSNP0), rep, nSNP0)
lower <- sapply(qnorm(maf, 0, 1),rep, nSNP0)
pha <- numeric(m)
hap.Gen <- matrix(2, nSNP0, nSNP0)
for(i in 1:nSNP0){
  upper[i, i] <- lower[i, i]
  lower[i, i] <- -Inf
  pha[i] <- mvtnorm::pmvnorm(lower=lower[i,], upper=upper[i, ], mean=rep(0, nSNP0), sigma=sigma)[1]
  hap.Gen[i, i] <- 1
}

pha1 <- mvtnorm::pmvnorm(lower=qnorm(maf, 0, 1), upper=Inf, mean=rep(0, nSNP0), sigma=sigma) 
hap.freq <- matrix(c(pha1[1], pha+(1-sum(pha)-pha1[1])/nSNP0), nrow=1)

hap <- paste(1, length(hap.freq), nSNP0, 1) # “1” for haplotypes, number of haplotypes, number of marker loci, number of times pattern repeated
hap.Gen <- rbind(rep(2, nSNP0), hap.Gen)
 
write.table(names, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(seed, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(traid, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(n.fam, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(n.snp, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(re, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(re.num, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(hap, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(hap.freq, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(hap.Gen, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
########################################


######################## To create the parameter file for genotype data #########
rho <- 0
nSNP0 <- 8
names <- paste("pedin_unrelated.pre",  paste0("pedfile_unrelated","_",rho, "_",nSNP0, ".pre")) ### name of pedigree file, name of output file
seed <- paste(23221, 1601, 21001)             ## three random seeds
traid <- 0   # number of columns for affection status/quantitative trait
n.fam <- 625000  # number of replicates/families
n.snp <- paste(nSNP0, 1)    # Total number of marker loci, # of times pattern to be repeated 
re <- 1  # “1” recomb fraction, “2” Kosambi map distance & “3” Haldane map distance  
re.num <- matrix(c(nSNP0-1, rep(0, nSNP0-1)), nrow=1)  # no recombination
G.num <- paste(2, nSNP0, 1)  # “2” genotype data will be provided, # of marker loci, # of times pattern repeated
maf <- runif(nSNP0, 0.001, 0.05)
G <- cbind(rep(2, nSNP0), maf)  ## Generate the genotype

dir <- "~/Rare-variant/GData"
write.table(names, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(seed, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(traid, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(n.fam, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(n.snp, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(re, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(re.num, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(G.num, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)
write.table(G, file=paste0(dir, "/input_unrelated","_", rho, "_",nSNP0,  ".dat"), quote=F, row.names = F, col.names = F, append=T)










