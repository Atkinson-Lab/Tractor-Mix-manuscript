args <- commandArgs(trailingOnly = TRUE)
No <- as.integer(args[1])


library(data.table)
SNP_meta = cbind(CHR = rep(1, 50000), POS = rep(".", 50000), SNP = rep(".", 50000), REF = rep("A", 50000), ALT = rep("G", 50000))


dirname = paste0("/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/Block_row", No, "/")
dir.create(dirname)
bLOCK_idx = ((No - 1) * 50000 + 1):(50000 * No)



GTot = fread("/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/GTot.tsv", sep = "\t")
GAfr = fread("/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/GAfr.tsv", sep = "\t")
GEur = fread("/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/GEur.tsv", sep = "\t")
LA = fread("/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/LA.tsv", sep = "\t")


write.table(cbind(SNP_meta, GTot[bLOCK_idx, ]),  paste0(dirname, "GTot.tsv"), sep = "\t", quote = F, row.names = F)
write.table(cbind(SNP_meta, GAfr[bLOCK_idx, ]),  paste0(dirname, "GAfr.tsv"), sep = "\t", quote = F, row.names = F)
write.table(cbind(SNP_meta, GEur[bLOCK_idx, ]),  paste0(dirname, "GEur.tsv"), sep = "\t", quote = F, row.names = F)
write.table(cbind(SNP_meta, LA[bLOCK_idx, ]),  paste0(dirname, "LA.tsv"), sep = "\t", quote = F, row.names = F)
