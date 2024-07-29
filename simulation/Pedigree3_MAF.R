source("Utils.R")
resdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/"
set.seed(123)
p = 1000000
MAFs = GenerateMafs(p, fst = 0.1, runif1 = 0.1, runif2 = 0.9)
write.table(cbind(SNP = paste0("SNP", 1:p), MAF_AFR =  MAFs$maf1, MAF_EUR = MAFs$maf2), 
	    paste0(resdir, "SNP_AF.tsv"),
	    quote = F, sep = "\t", col.names = F, row.names = F)
