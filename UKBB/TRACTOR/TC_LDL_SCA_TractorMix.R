args <- commandArgs(trailingOnly = TRUE)
CHR <- as.integer(args[1])

source("../Tractor-Mix/TractorMix.score.R")
infiles = c(paste0("~/Tractor-Mix/UKBB/data/TRACTOR/HAP_DOSE/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR", CHR, "-AC-phased.anc0.dosage.txt"),
            paste0("~/Tractor-Mix/UKBB/data/TRACTOR/HAP_DOSE/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR", CHR, "-AC-phased.anc1.dosage.txt"))

load("~/Tractor-Mix/UKBB/data/TRACTOR/NULLs/TC_Null.RData")
load("~/Tractor-Mix/UKBB/data/TRACTOR/NULLs/LDL_Null.RData")
load("~/Tractor-Mix/UKBB/data/TRACTOR/NULLs/SCA_Null.RData")

TractorMix.score(obj = SCA_Null,
                 infiles = infiles,
                 outfiles = paste0("~/Tractor-Mix/UKBB/res/SUMSTATS/SCA/Tractor-Mix1/sumstats_chr", CHR, ".tsv"), 
                 n_core = 8, chunk_size = 2048, AC_threshold = 50)

TractorMix.score(obj = TC_Null, 
                 infiles = infiles,
                 outfiles = paste0("~/Tractor-Mix/UKBB/res/SUMSTATS/TC/Tractor-Mix1/sumstats_chr", CHR, ".tsv"),
                 n_core = 8, chunk_size = 2048)


TractorMix.score(obj = LDL_Null, 
                 infiles = infiles,
                 outfiles = paste0("~/Tractor-Mix/UKBB/res/SUMSTATS/LDL/Tractor-Mix1/sumstats_chr", CHR, ".tsv"),
                 n_core = 8, chunk_size = 2048)





