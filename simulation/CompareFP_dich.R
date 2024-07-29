args <- commandArgs(trailingOnly = TRUE)
No <- as.integer(args[1])
paramNo <- as.integer(args[2])

library(GMMAT)
library(data.table)
library(MASS)
library(lmtest)
source("TractorMix.score.R")

datadir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/"
resdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/FP/Dich/"
restempdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/FP/Dich/temp/"
resblockdir = paste0(resdir, "Block_row", No, "/")

GTotdir = paste0(datadir, "Block_row", No, "/GTot.tsv")
GAfrdir = paste0(datadir, "Block_row", No, "/GAfr.tsv")
GEurdir = paste0(datadir, "Block_row", No, "/GEur.tsv")
LAdir = paste0(datadir, "Block_row", No, "/LA.tsv")
TractorMixinfilesdir = c(GAfrdir, GEurdir)

dir.create(resblockdir)

SampleInfodir = paste0(restempdir, "SampleInfo_param", paramNo, ".tsv")
model0_GMMATdir = paste0(restempdir, "model0_GMMAT_param", paramNo, ".Rdata")
model0_GENESISdir = paste0(restempdir, "model0_GENESIS_param", paramNo, ".Rdata")
model0_Truedir = paste0(restempdir, "model0_True_param", paramNo, ".Rdata")

GWAS_resdir = paste0(resblockdir, "GWAS_res_param", paramNo, ".tsv")
Tractor_resdir = paste0(resblockdir, "Tractor_res_param", paramNo, ".tsv")
GMMAT_resdir = paste0(resblockdir, "GMMAT_res_param", paramNo, ".tsv")
TractorMix_GMMAT_resdir = paste0(resblockdir, "TractorMix_GMMAT_res_param", paramNo, ".tsv")
GENESIS_resdir = paste0(resblockdir, "GENESIS_res_param", paramNo, ".tsv")
TractorMix_GENESIS_resdir = paste0(resblockdir, "TractorMix_GENESIS_res_param", paramNo, ".tsv")
True_resdir = paste0(resblockdir, "True_res_param", paramNo, ".tsv")
TractorMix_True_resdir = paste0(resblockdir, "TractorMix_True_res_param", paramNo, ".tsv")


SampleInfo = read.csv(SampleInfodir, sep = "\t")
load(model0_GMMATdir)
load(model0_GENESISdir)
load(model0_Truedir)


GTot = fread(GTotdir, sep = "\t")
GAfr = fread(GAfrdir, sep = "\t")
GEur = fread(GEurdir, sep = "\t")
LA = fread(LAdir, sep = "\t")

RunGWAS = function(row){
  geno = as.integer(row[6:length(row)])
  pval = summary(glm(SampleInfo$pheno ~ SampleInfo$PC1 + geno, family = binomial(link = "logit")))$coefficients[3, 4]
  return(c(row[1:5], pval = pval))
}


resGWAS = t(apply(GTot, 1, RunGWAS))
write.table(resGWAS, GWAS_resdir, quote = F, row.names = F, sep = "\t")

RunTractor = function(rowid){
  X = cbind(PC = SampleInfo$PC1, 
            LA = t(LA[rowid, 6:ncol(LA)]),
            GAfr = t(GAfr[rowid, 6:ncol(GAfr)]), 
            GEur = t(GEur[rowid, 6:ncol(GEur)]))
  
  m0 = glm(SampleInfo$pheno ~ X[,1:2], family = binomial(link = "logit"))
  m1 = glm(SampleInfo$pheno ~ X, family = binomial(link = "logit"))
  pval = lrtest(m0, m1)[2,5]
  return(c(GTot[rowid,1:5], pval = pval))
}

resTractor = t(sapply(1:nrow(GTot), RunTractor))
write.table(resTractor, Tractor_resdir, quote = F, row.names = F, sep = "\t")



############ GMMAT (GRM + PC) ###############
glmm.score(model0_GMMAT, 
             infile = GTotdir, 
             outfile = GMMAT_resdir, 
             infile.nrow.skip = 1,
             infile.ncol.skip = 5, 
             infile.ncol.print = 1:5,
             center = F,
             infile.header.print = c("CHR", "POS", "ID", "REF", "ALT"))

TractorMix.score(model0_GMMAT, TractorMixinfilesdir, TractorMix_GMMAT_resdir)



############ GENESIS (PCAir + PCRelate) ###############
glmm.score(model0_GENESIS, 
             infile = GTotdir, 
             outfile = GENESIS_resdir, 
             infile.nrow.skip = 1,
             infile.ncol.skip = 5, 
             infile.ncol.print = 1:5,
             center = F,
             infile.header.print = c("CHR", "POS", "ID", "REF", "ALT"))

TractorMix.score(model0_GENESIS, TractorMixinfilesdir, TractorMix_GENESIS_resdir)



############ True (Admprop + Kinship) ###############
glmm.score(model0_True, 
             infile = GTotdir, 
             outfile = True_resdir, 
             infile.nrow.skip = 1,
             infile.ncol.skip = 5, 
             infile.ncol.print = 1:5,
             center = F,
             infile.header.print = c("CHR", "POS", "ID", "REF", "ALT"))

TractorMix.score(model0_True, TractorMixinfilesdir, TractorMix_True_resdir)



