library(data.table)
Pheno_covariate_data <- fread(file = "~/PretermData/PLINK_filtered_dataset/PD_QC_MCpairs_phenotype.txt", header = TRUE, sep = " ")  ## 1626 samples have maternal information
BMINA_ID <- which(is.na(Pheno_covariate_data$maternal_bmi))
pheno_deleteNA <- Pheno_covariate_data[-BMINA_ID, ]

phenotype <- pheno_deleteNA$premie; table(phenotype)
BMI <- pheno_deleteNA$maternal_bmi; m_height <- pheno_deleteNA$maternal_height;  m_age <- pheno_deleteNA$maternal_age

load(file="~/PretermData/PLINK_filtered_dataset/PD_MGdata.RData")

PD_QC_MCpairs_pheno_BMI.bim <- read.table(file="~/PretermData/PLINK_filtered_dataset/PD_QC_MCpairs_pheno_BMI.bim")
all_SNP <- PD_QC_MCpairs_pheno_BMI.bim$V2
chro <- PD_QC_MCpairs_pheno_BMI.bim$V1
position <- PD_QC_MCpairs_pheno_BMI.bim$V4 

# Candidate gene or SNPs：RBPJ

#dist <- 300000 
#SNPs <- intersect(which(position >= (26105449+dist) & position <= (26435131+dist)) , which(chro==4)) 
#candi_SNPs <- all_SNP[SNPs]
#write.table(candi_SNPs, file = "~/PretermData/PLINK_filtered_dataset/RBPJ.txt", quote = FALSE, row.names = FALSE, col.names = FALSE,append=FALSE)

#ATP2A1, SEC16B
candi_SNPs <- as.vector(read.table(file = "~/PretermData/PLINK_filtered_dataset/SEC16B.txt", header = F))$V1
id <- which(all_SNP==candi_SNPs[1])
SNPs <- id:(id+length(candi_SNPs)-1)
  
# check if covariates NA
all(is.na(m_height) == FALSE) ## true
all(is.na(m_age) == FALSE) ## false so complete NA using mean values
m_age[is.na(m_age)] <- mean(m_age, na.rm = TRUE)

source("~/Rare-variant/Mixed_chi_squared_pvalue.R")
source("~/Rare-variant/aGE_functions.R")
source("~/Rare-variant/Statistic_functions1.R")
 
dir <- "~/Rare-variant/Real_data_result"

re <- list()
re$X <- cbind(BMI, m_height, m_age)
re$G <- as.matrix(MGdata)[, SNPs]
re$Y <- phenotype
p_inter <- Get_inter_glmm(re, alpha=0.1,pow = c(1:6, Inf), n.perm =1e4,
                                mainweights=poweights, interweights=poweights) 

p_joint <- Get_joint_glm(re, alpha=0.1,pow = c(1:6, Inf), n.perm =1e4,
                            mainweights=poweights, interweights=poweights) 

p.aGE.inter <- p_inter$p.aGE
p.IndaGE.inter <- p_inter$p.IndaGE
p.waGE.inter <- p_inter$p.waGE
p.rareGE.inter <- p_inter$p.rareGE
p.IndrareGE.inter <- p_inter$p.IndrareGE
p.wrareGE.inter <- p_inter$p.wrareGE
p0.inter <- p_inter$p0
#######################################
p.aGE.joint <- p_joint$p.aGE
p.IndaGE.joint <- p_joint$p.IndaGE
p.waGE.joint <- p_joint$p.waGE
p.rareGE.joint <- p_joint$p.rareGE
p.IndrareGE.joint <- p_joint$p.IndrareGE
p.wrareGE.joint <- p_joint$p.wrareGE
p0.joint <- p_joint$p0

pval_inter <- rbind(cbind(p.aGE.inter, p.IndaGE.inter, p.waGE.inter), cbind(p.rareGE.inter, p.IndrareGE.inter, p.wrareGE.inter))
pval_joint <- rbind(cbind(p.aGE.joint, p.IndaGE.joint, p.waGE.joint), cbind(p.rareGE.joint, p.IndrareGE.joint, p.wrareGE.joint))

write.table(pval_inter, file = paste0(dir, "/inter_result_SEC16B.txt"), quote = FALSE,append=FALSE)
write.table(pval_joint, file = paste0(dir, "/joint_result_SEC16B.txt"), quote = FALSE,append=FALSE)
write.table(c(p0.inter, p0.joint), file = paste0(dir, "/ind_SEC16B.txt"), quote = FALSE,append=FALSE)

