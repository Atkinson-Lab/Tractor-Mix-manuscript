args <- commandArgs(trailingOnly = TRUE)
No <- as.integer(args[1])

# This script is used to simulate statistical power of TractorMix model for continuous phenotypes
# We also use this script to evaluate the effect sizes estimates from the wald test and score test
# We will test 24 parameters, with effect size changes (AFR or/and EUR)
# The simulation uses true pedigree and admixture proportion instead of output from GENESIS
# In realistic setting, the performance is likely less good
# We continue to use pedigree3 as the true kinship
# The script will produce a tsv file for each row, with p values and effect size estimates recorded
# The results can be summerize into power curve and effect sizes estimates

source("Utils.R")
source("Pedigree3.R")
source("TractorMix.score.R")
source("TractorMix.wald.R")
source("Utils.R")
source("Pedigree3.R")
library(kinship2)
library(GMMAT)
library(MASS)

resdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/Power_cont/"
dir.create(paste0(resdir, "data/"))
resdir_data = paste0(resdir, "data/", paste0("data", No), "/")
resdir_res = paste0(resdir, "res/")
resdir_plot = paste0(resdir, "plot/")


dir.create(resdir)
dir.create(resdir_data)
dir.create(resdir_res)
dir.create(resdir_plot)

infilesTractorMix = c(paste0(resdir_data, "AFRg.txt"), paste0(resdir_data, "EURg.txt"))
outTractorMix = paste0(resdir_data, "TractorMix.res")
infilesGMMAT = paste0(resdir_data, "TOTg.txt")

##############################################################################################
# Tractor mix was run n times, each time with fixed parameters but with different phenotypes #
##############################################################################################
RunTractorMix <- function(i, n, infiles, outfile, X, betas, tau2, Kin, sigma2, Admprop){
  # occasionally results can contain NAs, and this can make the dimention incorrect
  res = rep(NA, 6)
  while (anyNA(res)){
      pheno = mvrnorm(1, X %*% betas, tau2 * Kin + sigma2 * diag(n))
      PHENO = data.frame(Admprop, pheno, ID = paste0("Sample", 1:n))
      model0_gmmat_true = glmmkin(fixed = pheno ~  1 + Admprop, data = PHENO, id = "ID", kins = Kin, family = gaussian(link = "identity"))
      TractorMix.score(model0_gmmat_true, infiles, outfile)
      TractorMix_res = read.csv(outfile, sep = "\t")
      
      Tractormix.wald = TractorMix.wald(infiles = infiles, 
                                 data = PHENO, 
                                 fixed = pheno~ Admprop, kin = Kin, id = ID, family = gaussian(link = "identity"))
    
      
      glmm.score(model0_gmmat_true, 
                 infile = paste0(resdir_data, "TOTg.txt"), 
                 outfile = paste0(resdir_data, "TOTg.res"), 
                 infile.nrow.skip = 1,
                 infile.ncol.skip = 5, 
                 infile.ncol.print = 1:5,
                 center = F,
                 infile.header.print = c("CHR", "POS", "ID", "REF", "ALT"))
      
      # once a while GMMAT will fail, handing exception
      gmmatP = try(read.csv(paste0(resdir_data, "TOTg.res"),sep = "\t")$PVAL)
      if (class(file) == "try-error") {
            cat("Caught an error during read gmmat res")
            gmmatP <- NA
      }
      res = c(gmmatP, 
               TractorMix_res$P, 
               TractorMix_res$Eff_anc1, 
               TractorMix_res$Eff_anc2, 
               Tractormix.wald$AFRg$BETA, 
               Tractormix.wald$EURg$BETA)
  }
  
  return(res)
}






# params is a vector of parameters
# params[1] is the sample size
# params[2:3] is the distribution of global admixture 
# params[4] is the allele frequency for AFR
# params[5] is the allele frequency for EUR
# params[6:9] is betas (intercept, beta_adm, beta_AFR, beta_EUR)
# params[10] is sigma2
# params[11] is tau2
# params[12] is the number of iterations 
# params[13] is the number of parameter space 

# The kinship is set up in a way that half of individuals are independnet and half are relatives

RunSimulation <- function(params){

  n = params[1]
  if ( (n/2) %% 10 != 0  ) stop('Genotype seperation error')
  Admprop = rbeta(params[1], params[2], params[3])
  AFRmaf = params[4]
  EURmaf = params[5]
  AFRla = sapply(Admprop, function(admprop){rbinom(1, size = 2, prob = admprop)} )
  EURla = 2 - AFRla
  AFRg = matrix(sapply(AFRla, function(la){rbinom(1, size = la, prob = AFRmaf)}))
  EURg = matrix(sapply(EURla, function(la){rbinom(1, size = la, prob = EURmaf)}))
  TOTg = AFRg + EURg
  
  WriteGenotype(AFRg, paste0(resdir_data, "AFRg.txt"))
  WriteGenotype(EURg, paste0(resdir_data, "EURg.txt"))
  WriteGenotype(TOTg, paste0(resdir_data, "TOTg.txt"))
  
  X = cbind(1, Admprop, AFRg, EURg)
  betas = params[6:9]
  sigma2 = params[10]
  tau2 = params[11]
  Kin = MakeGRM(params[1]/2, params[1]/2)
  
  res = matrix(unlist(lapply(1:params[12], RunTractorMix, n, infilesTractorMix, outTractorMix, X, betas, tau2, Kin, sigma2, Admprop)), ncol = 6, byrow = T)
  colnames(res) = c("gmmat_P", "Tractor_score_P", 
                    "Tractor_score_AFR_beta", "Tractor_score_EUR_beta",
                    "Tractor_wald_AFR_beta", "Tractor_wald_EUR_beta")
  
  write.table(x = res, file = paste0(resdir_res, "res", params[13], ".tsv" ), quote = F, col.names = T, row.names = F, sep = "\t" )
}


paramsdf = read.csv("Power_params_cont.tsv", sep = "\t")[No,]

######## run program ####################
apply(paramsdf, 1, RunSimulation)
