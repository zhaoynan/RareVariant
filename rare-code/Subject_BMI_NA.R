# keep subjects whose BMI exist

library(data.table)
Pheno_covariate_data <- fread(file = "~/PretermData/PLINK_filtered_dataset/PD_QC_MCpairs_phenotype.txt", header = TRUE, sep = " ")  ## 1626 samples have maternal information
PD_QC.fam <- read.table(file="~/PretermData/PLINK_QC_process/PD_QC.fam")
BMINA_ID <- which(is.na(Pheno_covariate_data$maternal_bmi))
pheno_deleteNA <- Pheno_covariate_data[-BMINA_ID, ]

phenotype <- pheno_deleteNA$premie; table(phenotype)
BMI <- pheno_deleteNA$maternal_bmi; m_height <- pheno_deleteNA$maternal_height;  m_age <- pheno_deleteNA$maternal_age
m_ID <- pheno_deleteNA$m_genevaID; c_ID <- pheno_deleteNA$c_genevaID

subject.keep <- cbind(c(PD_QC.fam$V1[match(m_ID, PD_QC.fam$V2)], PD_QC.fam$V1[match(c_ID, PD_QC.fam$V2)]), c(m_ID, c_ID)) 
write.table(subject.keep, file = "~/PretermData/PLINK_filtered_dataset/subject.keep.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ./plink --bfile ~/PretermData/PLINK_QC_process/PD_QC --keep ~/PretermData/PLINK_filtered_dataset/subject.keep.txt --make-bed --out ~/PretermData/PLINK_filtered_dataset/PD_QC_MCpairs_pheno_BMI
#----------------------------------------------------------------------------------------#

#  .bed and .fam data were converted into 0,1,2 raw data
#--------------------------------PLINK---------------------------------------------------#
# ./plink --bfile ~/PretermData/PLINK_filtered_dataset/PD_QC_MCpairs_pheno_BMI --recodeA --out ~/PretermData/PLINK_filtered_dataset/PD_QC_MCpairs_additive
#----------------------------------------------------------------------------------------#
