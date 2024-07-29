library(data.table)
file1 = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/Block_row1/GTot.tsv"
file2 = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/Block_row2/GTot.tsv"
GRMdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/GRM.tsv"
PCdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/PC.tsv"

GTot_blk1 = fread(file1, sep = "\t")
GTot_blk2 = fread(file2, sep = "\t")

mat = t(rbind(GTot_blk1[,6:ncol(GTot_blk1)], GTot_blk2[,6:ncol(GTot_blk2)]))
af = colMeans(mat)/2
mat_scale = sapply(1:ncol(mat), function(i){(mat[,i] - 2 * af[i])/sqrt(2 * af[i] * (1 - af[i]))} )
GRM = mat_scale %*% t(mat_scale)/(ncol(mat))
PCs = prcomp(mat, scale = T)$x[,1:5]


write.table(GRM, GRMdir,row.names = F, col.names=F, sep = "\t")
write.table(PCs, PCdir, quote = F, sep = "\t", row.names = F)
