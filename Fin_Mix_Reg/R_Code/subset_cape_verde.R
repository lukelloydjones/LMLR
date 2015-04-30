# ==============================================================================
# Script to annotate the SNPs found by Bayes R to help make biological sense
# of the positions found
# ==============================================================================
source("~/Dropbox/Post_Doc_QBI/R_Functions/plink_read_func.R")
geno <- read.plink("~/Desktop/Cape_Verde_BayesR/Original_Raw_Data/CapeVerdeFinal")
bim  <- read.table("~/Desktop/Cape_Verde_BayesR/Original_Raw_Data/CapeVerdeFinal.bim", fill = T)
fam  <- read.table("~/Desktop/Cape_Verde_BayesR/Original_Raw_Data/CapeVerdeFinal.fam", header = F)
rownames(geno) <- fam[, 2]
colnames(geno) <- bim[, 2]
geno.2 <- geno[, sample(seq(1, 903837), 10000)]
geno.2[1:10, 1:10]
write.table(geno.2, "cape_verde_r10k.txt", row.names = T, col.names = T, sep = "\t", quote = F)
pheno <- read.table("~/Desktop/Cape_Verde_BayesR/Original_Raw_Data/CapeVerdeFinal_Pheno.fam", header = F)
plot(density(pheno[, 7]))
colnames(pheno) <- c("FID", "IID", "P1", "P2" ,"P3" ,"EYE" ,"SKIN", "P4" , "P5" )
head(pheno)
pheno.2 <- pheno[, -c(1,3,4,5, 8, 9)]
head(pheno.2)
write.table(pheno.2, "cape_verde_pheno.txt", row.names = F, col.names = T, sep = "\t", quote = F)