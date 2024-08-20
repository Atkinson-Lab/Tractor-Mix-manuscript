args <- commandArgs(trailingOnly = TRUE)
CHR <- as.integer(args[1])

source("../Tractor-Mix/TractorMix.score.update.R")
infiles = c(paste0("~/Tractor-Mix/UKBB/data/TRACTOR/HAP_DOSE/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR", CHR, "-AC-phased.anc0.dosage.txt"),
            paste0("~/Tractor-Mix/UKBB/data/TRACTOR/HAP_DOSE/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR", CHR, "-AC-phased.anc1.dosage.txt"))

load("~/Tractor-Mix/UKBB/data/TRACTOR/NULLs/SCA_Null.RData")

TractorMix.score(obj = SCA_Null,
                 infiles = infiles,
                 outfiles = paste0("~/Tractor-Mix/UKBB/res/SUMSTATS/SCA/Tractor-Mix/sumstats_chr", CHR, "_update.tsv"), 
                 n_core = 8, chunk_size = 2048)





