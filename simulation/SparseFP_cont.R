args <- commandArgs(trailingOnly = TRUE)
No <- as.integer(args[1])

library(kinship2)
library(GMMAT)
library(MASS)
source("TractorMix.sparse.score.R")
source("Pedigree3.R")
FullKinship = MakeGRM(2500, 2500)
SparseKinship = Matrix(FullKinship, sparse = TRUE)



data = c(paste0("/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/Block_row", No, "/GAfr.tsv"), 
         paste0("/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/Block_row", No, "/GEur.tsv"))

admdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Pedigree3/data/"
typedir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Sparse/Cont/"
blockdir = paste0(typedir, "Block_row", No, "/")
dir.create(typedir)
dir.create(blockdir)



Admprop = read.csv(paste0(admdir, "Admprop.tsv"), sep = "\t", header = F)
colnames(Admprop) = c("ID", "Admprop")
beta0 = 1
beta1 = 3
tau2 = 3
sigma2 = 1
niter = 2000

for (phenoNo in 1:niter){
  resdir = paste0(blockdir, "res", phenoNo, ".tsv")
  
  SampleInfo = cbind(Admprop, pheno = mvrnorm(1, 1 * beta0 + Admprop[,2] * beta1,  tau2 * FullKinship + sigma2 * diag(5000)))
  model0_SpaseKinship =  glmmkin(pheno ~ Admprop, data = SampleInfo, 
                                 kin = SparseKinship, id = "ID", family = gaussian(link = "identity"))
  
  TractorMix.score(model0_SpaseKinship, data, resdir)
}

