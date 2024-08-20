library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(data.table)
file_prefix = "~/Tractor-Mix/UKBB/data/TRACTOR/GENESIS/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-prune"

snpgdsBED2GDS(bed.fn = paste0(file_prefix, ".bed"), 
              bim.fn = paste0(file_prefix, ".bim"), 
              fam.fn = paste0(file_prefix, ".fam"), 
              out.gdsfn = paste0(file_prefix, ".gds"))

KING <- fread(paste0(file_prefix, ".king"), data.table = FALSE)
KING.ID <- fread(paste0(file_prefix, ".king.id"), data.table = FALSE)
KINGmat = as.matrix(KING)
row.names(KINGmat) = KING.ID[,1]
colnames(KINGmat) = KING.ID[,1]  

UKBB_geno <- GdsGenotypeReader(filename = paste0(file_prefix, ".gds"))
UKBB_genoData <- GenotypeData(UKBB_geno)
pcair_res <- pcair(UKBB_genoData, kinobj = KINGmat, divobj = KINGmat, 
                   num.cores = 16)



### PC-Relate
UKBB_genoData <- GenotypeBlockIterator(UKBB_genoData)
pcrelate_res <- pcrelate(UKBB_genoData, pcs = pcair_res$vectors[,1:2], 
                         training.set = pcair_res$unrels,
                         BPPARAM = BiocParallel::SerialParam(),
                         ibd.probs = FALSE)


iids <- as.character(getScanID(UKBB_genoData))

# process PCAirmat, a matrix with 10 columns 
PCAirmat = pcair_res$vectors[iids,1:5]
colnames(PCAirmat) = paste0("PC", 1:5)

# process PCRelatemat, a sparse matrix 
PCRelatemat = pcrelateToMatrix(pcrelate_res, sample.include = iids, scaleKin = 2)[iids,iids]
PCRelatemat_sparse = PCRelatemat
PCRelatemat_sparse[PCRelatemat_sparse < 0.05] = 0
PCRelatemat_sparse = as(PCRelatemat_sparse, "sparseMatrix") 

### check if id matches
all.equal(row.names(PCAirmat), row.names(PCRelatemat_sparse), colnames(PCRelatemat_sparse))

PC_GRMs = list(PCAirmat = PCAirmat, 
               PCRelatemat = PCRelatemat_sparse)

save(PC_GRMs, file = paste0(file_prefix, "-PC_GRMs.RData"))
