args <- commandArgs(trailingOnly = TRUE)
CHR <- as.integer(args[1])

library(GMMAT)
load("~/Tractor-Mix/UKBB/data/TRACTOR/NULLs/TC_Null.RData")
load("~/Tractor-Mix/UKBB/data/TRACTOR/NULLs/LDL_Null.RData")
load("~/Tractor-Mix/UKBB/data/TRACTOR/NULLs/SCA_Null.RData")

geno.file = paste0("UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-chr", CHR)


glmm.score(TC_Null, infile = geno.file, center = F, missing.method = "omit",
           outfile = paste0("TC_sumstats_chr", CHR, ".tsv"))

glmm.score(SCA_Null, infile = geno.file, center = F,missing.method = "omit",
           outfile = paste0("SCA_sumstats_chr", CHR, ".tsv"))

glmm.score(LDL_Null, infile = geno.file, center = F,missing.method = "omit",
           outfile = paste0("LDL_sumstats_chr", CHR, ".tsv"))

