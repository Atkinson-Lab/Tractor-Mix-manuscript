args <- commandArgs(trailingOnly = TRUE)
No <- as.integer(args[1])

source("Utils.R")
resdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/"
Blockdir = paste0(resdir, "Block", No, "/")
dir.create(Blockdir)
MAFs = read.csv("/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/SNP_AF.tsv", sep = "\t", header = F)
n = 100
SampleID = paste0("Sample", ((No - 1) * n + 1):((No - 1) * n + n))


GTot_Ind = matrix(data = NA, nrow = n, ncol = nrow(MAFs))
GAfr_Ind = matrix(data = NA, nrow = n, ncol = nrow(MAFs))
GEur_Ind = matrix(data = NA, nrow = n, ncol = nrow(MAFs))
LA_Ind = matrix(data = NA, nrow = n, ncol = nrow(MAFs))
Admprop_Ind = matrix(data = NA, nrow = n, ncol = 1)

row.names(GTot_Ind) = SampleID
row.names(GAfr_Ind) = SampleID
row.names(GEur_Ind) = SampleID
row.names(LA_Ind) = SampleID
row.names(Admprop_Ind) = SampleID


for (i in 1:n){
    param1 = sample(c(8,12), size = 1)
    ind = GenerateAdm(MAFs[,2], MAFs[,3],  rbeta(1, param1, 20 - param1))
    GTot_Ind[i,] = ind$hap1 + ind$hap2
    GAfr_Ind[i,] = GetAncestrySpecCount(ind, 0)
    GEur_Ind[i,] = GetAncestrySpecCount(ind, 1)
    LA_Ind[i,]   = ind$la1 + ind$la2
    Admprop_Ind[i, 1] = GetGlobAncestry(ind)
}


write.table(t(GTot_Ind), paste0(Blockdir, "GTot.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(t(GAfr_Ind), paste0(Blockdir, "GAfr.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(t(GEur_Ind), paste0(Blockdir, "GEur.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(t(LA_Ind), paste0(Blockdir, "LA.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(Admprop_Ind, paste0(Blockdir, "Admprop.tsv"), quote = F, sep = "\t", col.names = F, row.names = T)


