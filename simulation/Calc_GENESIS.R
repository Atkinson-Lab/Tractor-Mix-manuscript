library(data.table)
library(genio)
library(GENESIS)
library(SNPRelate)
library(GWASTools)

file1 = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/Block_row1/GTot.tsv"
file2 = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/Block_row2/GTot.tsv"
tempdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/temp/GTot"
PCAirdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/PCAir.tsv"
PCRelatedir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/PCRelate.tsv"

dir.create("/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/temp/")


GTot_blk1 = fread(file1, sep = "\t")
GTot_blk2 = fread(file2, sep = "\t")

mat = t(rbind(GTot_blk1[,6:ncol(GTot_blk1)], GTot_blk2[,6:ncol(GTot_blk2)]))

bim = as.data.frame(cbind(chr = rep(1, ncol(mat)), id = ".", posg = 1:ncol(mat), pos = 1:ncol(mat), ref = "C", alt = "G"))
fam = as.data.frame(cbind(fam = paste0("Sample", 1:nrow(mat)), id = paste0("Sample", 1:nrow(mat)), pat = ".", mat = ".", sex = ".", pheno = "." ))
write_plink(tempdir, t(mat), bim = bim, fam = fam)

snpgdsBED2GDS(bed.fn = paste0(tempdir, ".bed"), 
              bim.fn = paste0(tempdir, ".bim"), 
              fam.fn = paste0(tempdir, ".fam"), 
              out.gdsfn = paste0(tempdir, ".gds"))

gds <- snpgdsOpen(paste0(tempdir, ".gds"))
KINGmat = snpgdsIBDKING(gds)$kinship
snpgdsClose(gds)

gds <- GdsGenotypeReader(filename = paste0(tempdir, ".gds"))
gdsData <- GenotypeData(gds)
row.names(KINGmat) = paste0("Sample", 1:nrow(KINGmat))
colnames(KINGmat) = paste0("Sample", 1:ncol(KINGmat))


# PC Air calculation 
PCair <- pcair(gdsData, kinobj = KINGmat, divobj = KINGmat)
PCairpart <- pcairPartition(kinobj = KINGmat, divobj = KINGmat)

# PC relate calculation 
gdsData <- GenotypeBlockIterator(gdsData)
PCrelate <- pcrelate(gdsData, pcs = PCair$vectors[,1:2], 
                          training.set = PCairpart$unrels,
                          BPPARAM = BiocParallel::SerialParam())

GRM_PC = pcrelateToMatrix(PCrelate, sample.include = paste0("Sample", 1:nrow(KINGmat)), scaleKin = 2)
PCrelate_out = as.matrix(GRM_PC)[paste0("Sample", 1:nrow(KINGmat)), paste0("Sample", 1:ncol(KINGmat))]

PCair_out = PCair$vectors[,1:5]
colnames(PCair_out) = paste0("PC", 1:5)

write.table(PCair_out, PCAirdir, quote = F, sep = "\t", row.names = F)
write.table(PCrelate_out, PCRelatedir, row.names = F, col.names=F, sep = "\t")
