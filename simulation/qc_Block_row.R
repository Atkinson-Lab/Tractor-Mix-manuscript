args <- commandArgs(trailingOnly = TRUE)
No <- as.integer(args[1])

library(data.table)
dirname = paste0("/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/Block_row", No, "/")


Filter = function(row){
  af = mean(row)/2
  return(af > 0.05 &  af < 0.95)
}


GTot = fread(paste0(dirname, "GTot.tsv"), sep = "\t")
GAfr = fread(paste0(dirname, "GAfr.tsv"), sep = "\t")
GEur = fread(paste0(dirname, "GEur.tsv"), sep = "\t")
LA = fread(paste0(dirname, "LA.tsv"), sep = "\t")

# Use Allele frequency as a filter
If_keep  = apply(GTot[,6:ncol(GTot)], 1, Filter)

write.table(GTot[If_keep, ],  paste0(dirname, "GTot.tsv"), sep = "\t", quote = F, row.names = F)
write.table(GAfr[If_keep, ],  paste0(dirname, "GAfr.tsv"), sep = "\t", quote = F, row.names = F)
write.table(GEur[If_keep, ],  paste0(dirname, "GEur.tsv"), sep = "\t", quote = F, row.names = F)
write.table(LA[If_keep, ],  paste0(dirname, "LA.tsv"), sep = "\t", quote = F, row.names = F)
