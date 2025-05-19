library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(dplyr)
library(data.table)
library(zoo)
library(cli)
library(Matrix)
library(withr)
library(BiocParallel)

# Define custom paste operator
"%&%" <- function(a, b) paste0(a, b)

# Set up project directories
proj_dir <- "/well/emberson/users/pom143/projects/LAI_TractorMix/"
out_dir <- proj_dir %&% "output-files/"

# Define paths to data files in out_dir
gds_path <- out_dir %&% "analysis-set-1_bmi.gds"
snpset_path <- out_dir %&% "analysis-set-1_bmi_snpset.RData"
pcair_results_path <- out_dir %&% "analysis-set-1_bmi_pcair_results.RData"

# Check if files exist
if (!file.exists(gds_path)) stop("GDS file not found in out_dir.")
if (!file.exists(snpset_path)) stop("SNP set file not found in out_dir.")
if (!file.exists(pcair_results_path)) stop("PCA results file not found in out_dir.")

# Set SNP block size and parallel parameters
snpBlockSize <- 5000  # SNP block size
bpparam <- BiocParallel::SnowParam(workers = 4, type = "SOCK")

# Load SNP set and PCA results
load(snpset_path)  # loads snpset
load(pcair_results_path)  # loads pcair_res

pruned <- unlist(snpset, use.names=FALSE)
total_snps <- length(pruned)

# Load GDS data with SNP block size
geno <- GdsGenotypeReader(gds_path)
genoData <- GenotypeData(geno)
genoData <- GenotypeBlockIterator(genoData, snpInclude=pruned, snpBlock=snpBlockSize)

# Define unrelated samples and PCs
unrel_set <- pcair_res$unrels
pcs <- pcair_res$vectors[, 1:2]

# Run PCRelate on the entire unrel_set without sample selection
pcrelate_res2 <- pcrelate(genoData, pcs = pcs, training.set = unrel_set, verbose = TRUE, ibd.probs = FALSE, 
                     BPPARAM = bpparam)
                     

save(pcrelate_res2, file = out_dir %&% "analysis-set-1_bmi_pcrelate_results_33500.RData")

# Close the GDS file after the analysis
close(geno)

