args <- commandArgs(trailingOnly = TRUE)
No <- as.integer(args[1])

library(GMMAT)
library(data.table)
library(MASS)

datadir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/"
FPdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/FP/"
resdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/FP/Dich/"
restempdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/FP/Dich/temp/"

dir.create(FPdir)
dir.create(resdir)
dir.create(restempdir)

SampleInfodir = paste0(restempdir, "SampleInfo_param", No, ".tsv")
model0_GMMATdir = paste0(restempdir, "model0_GMMAT_param", No, ".Rdata")
model0_GENESISdir = paste0(restempdir, "model0_GENESIS_param", No, ".Rdata")
model0_Truedir = paste0(restempdir, "model0_True_param", No, ".Rdata")


# the beta0 is determined based on prevalence constraints (ignore tau2)
params = as.numeric(read.csv("FP_params_dich.tsv", sep = "\t")[No,])

FitNull <- function(params){
  beta0 = params[1]
  beta1 = params[2]
  tau2 = params[3]
  prevalence = params[4]
  
  
  Admprop = read.csv(paste0(datadir, "Admprop.tsv"), sep = "\t", header = F)
  colnames(Admprop) = c("ID", "Admprop")
  PC = read.csv(paste0(datadir, "PC.tsv"), sep = "\t", header = T)
  PCAir = read.csv(paste0(datadir, "PCAir.tsv"), sep = "\t", header = T)
  colnames(PCAir) = paste0("PCAir", 1:5)
  SampleInfo = cbind(Admprop, PC, PCAir)
  
  Kinship = as.matrix(fread(paste0(datadir, "Kinship.tsv"), sep = "\t", header = F))
  GRM = as.matrix(fread(paste0(datadir, "GRM.tsv"), sep = "\t", header = F))
  PCRelate = as.matrix(fread(paste0(datadir, "PCRelate.tsv"), sep = "\t", header = F))
  
  row.names(Kinship) = paste0("Sample", 1:5000)
  colnames(Kinship) = paste0("Sample", 1:5000)
  row.names(GRM) = paste0("Sample", 1:5000)
  colnames(GRM) = paste0("Sample", 1:5000)
  row.names(PCRelate) = paste0("Sample", 1:5000)
  colnames(PCRelate) = paste0("Sample", 1:5000)
  
  probs = 1/(1 + exp(-mvrnorm(1, 1 * beta0 + SampleInfo$Admprop * beta1, tau2 * Kinship)))
  SampleInfo$pheno = sapply(probs, function(prob){rbinom(n = 1, size = 1, prob = prob)} )

  write.table(SampleInfo, SampleInfodir, sep = "\t", quote = F, row.names = F)
  
  
  model0_GMMAT =   glmmkin(fixed = pheno ~  1 + PC1, data = SampleInfo, id = "ID", kins = GRM, family = binomial(link = "logit"))
  model0_GENESIS = glmmkin(fixed = pheno ~  1 + PCAir1, data = SampleInfo, id = "ID", kins = PCRelate, family = binomial(link = "logit"))
  model0_True =    glmmkin(fixed = pheno ~  1 + Admprop, data = SampleInfo, id = "ID", kins = Kinship, family = binomial(link = "logit"))
  
  save(model0_GMMAT,file = model0_GMMATdir)
  save(model0_GENESIS,file = model0_GENESISdir)
  save(model0_True, file = model0_Truedir)
}

FitNull(params)


