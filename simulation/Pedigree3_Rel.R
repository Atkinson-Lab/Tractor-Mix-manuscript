args <- commandArgs(trailingOnly = TRUE)
No <- as.integer(args[1])

source("Utils.R")
source("Pedigree3.R")
resdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/"
Blockdir = paste0(resdir, "Block", No, "/")
dir.create(Blockdir)
MAFs = read.csv("/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/SNP_AF.tsv", sep = "\t", header = F)
n = 100
SampleID = paste0("Sample", ((No - 1) * n + 1):((No - 1) * n + n))


GTot_Rel = matrix(data = NA, nrow = n, ncol = nrow(MAFs))
GAfr_Rel = matrix(data = NA, nrow = n, ncol = nrow(MAFs))
GEur_Rel = matrix(data = NA, nrow = n, ncol = nrow(MAFs))
LA_Rel = matrix(data = NA, nrow = n, ncol = nrow(MAFs))
Admprop_Rel = matrix(data = NA, nrow = n, ncol = 1)

row.names(GTot_Rel) = SampleID
row.names(GAfr_Rel) = SampleID
row.names(GEur_Rel) = SampleID
row.names(LA_Rel) = SampleID
row.names(Admprop_Rel) = SampleID


MAFs = list(maf1 = MAFs[,2], 
            maf2 = MAFs[,3])
famsize = 10 

for (i in 1:(n/famsize)){
  param1 = sample(c(8,12), size = 1)
  fam = MakePedigree(MAFs, param1)
  GTot_Rel[(famsize * (i - 1) + 1):(famsize * i),] = fam$GenoMatTot
  GAfr_Rel[(famsize * (i - 1) + 1):(famsize * i),] = fam$GenoMatAfr
  GEur_Rel[(famsize * (i - 1) + 1):(famsize * i),] = fam$GenoMatEur
  LA_Rel[(famsize * (i - 1) + 1):(famsize * i),]   = fam$LAMat
  Admprop_Rel[(famsize * (i - 1) + 1):(famsize * i),1] = fam$Admprop
}


write.table(t(GTot_Rel), paste0(Blockdir, "GTot.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(t(GAfr_Rel), paste0(Blockdir, "GAfr.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(t(GEur_Rel), paste0(Blockdir, "GEur.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(t(LA_Rel), paste0(Blockdir, "LA.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(Admprop_Rel, paste0(Blockdir, "Admprop.tsv"), quote = F, sep = "\t", col.names = F, row.names = T)

