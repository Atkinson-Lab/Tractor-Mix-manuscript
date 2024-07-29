source("Utils.R")
source("Pedigree3.R")
source("TractorMix.ancP.score.R")
source("TractorMix.wald.R")
library(kinship2)
library(GMMAT)
library(MASS)

resdir = "/storage/atkinson/home/u241265/Tractor-Mix/simu/ancP/"
dir.create(resdir)



n = 1000
betaparam1 = 8
betaparam2 = 2
AFRmaf = 0.1
EURmaf = 0.2
betas_null_cont = c(1,2,0,0)
betas_alt_cont = c(1,2,1,0.5)
sigma2 = 1
tau2 = 1
niter = 1000
Kin = MakeGRM(n/2, n/2)


############ Continuous #####################
RunSimulation_cont <- function(No, betas){
  resdir_data = paste0(resdir, paste0("data", No), "/")
  dir.create(resdir_data)
  infilesTractorMix = c(paste0(resdir_data, "AFRg.txt"), paste0(resdir_data, "EURg.txt"))
  outTractorMix = paste0(resdir_data, "TractorMix.res")
  
  Admprop = rbeta(n, betaparam1, betaparam2)
  AFRla = sapply(Admprop, function(admprop){rbinom(1, size = 2, prob = admprop)} )
  EURla = 2 - AFRla
  AFRg = matrix(sapply(AFRla, function(la){rbinom(1, size = la, prob = AFRmaf)}))
  EURg = matrix(sapply(EURla, function(la){rbinom(1, size = la, prob = EURmaf)}))
  
  WriteGenotype(AFRg, paste0(resdir_data, "AFRg.txt"))
  WriteGenotype(EURg, paste0(resdir_data, "EURg.txt"))
  
  X = cbind(1, Admprop, AFRg, EURg)
  
  res = rep(NA, 4)
  while (anyNA(res)){
    pheno = mvrnorm(1, X %*% betas, tau2 * Kin + sigma2 * diag(n))
    PHENO = data.frame(Admprop, pheno, ID = paste0("Sample", 1:n))
    
    model0_gmmat_true = glmmkin(fixed = pheno ~  1 + Admprop, data = PHENO, id = "ID", kins = Kin, family = gaussian(link = "identity"))
    TractorMix.score(model0_gmmat_true, infilesTractorMix, outTractorMix)
    TractorMix_score_res = as.numeric(read.csv(outTractorMix, sep = "\t")[,10:11])
    
    Tractormix.wald_restbl = TractorMix.wald(infiles = infilesTractorMix, 
                                     data = PHENO, 
                                     fixed = pheno~ Admprop, kin = Kin, id = ID, family = gaussian(link = "identity"))
    
    TractorMix_wald_res = c(Tractormix.wald_restbl$AFRg[,"P"], 
                            Tractormix.wald_restbl$EURg[,"P"])
    
    res = c(TractorMix_score_res, TractorMix_wald_res)
    names(res) = c("AFR_score", "EUR_score", "AFR_wald", "EUR_wald")
  }
  file.remove(paste0(resdir_data, "AFRg.txt"), paste0(resdir_data, "EURg.txt"), paste0(resdir_data, "TractorMix.res"))
  return(res)
}




pval_res_null_cont =  t(sapply(1:1000, RunSimulation_cont, betas_null_cont))
pval_res_alt_cont =  t(sapply(1:1000, RunSimulation_cont, betas_alt_cont))




############ Dichotomous #####################
beta0space = seq(-4, 0, 0.1)
prevalence = 0.2
betas_null_dich = c(2,0,0)
betas_alt_dich = c(2,2,1)

RunSimulation_dich <- function(No, betas){
  resdir_data = paste0(resdir, paste0("data", No), "/")
  dir.create(resdir_data)
  infilesTractorMix = c(paste0(resdir_data, "AFRg.txt"), paste0(resdir_data, "EURg.txt"))
  outTractorMix = paste0(resdir_data, "TractorMix.res")
  
  Admprop = rbeta(n, betaparam1, betaparam2)
  AFRla = sapply(Admprop, function(admprop){rbinom(1, size = 2, prob = admprop)} )
  EURla = 2 - AFRla
  AFRg = matrix(sapply(AFRla, function(la){rbinom(1, size = la, prob = AFRmaf)}))
  EURg = matrix(sapply(EURla, function(la){rbinom(1, size = la, prob = EURmaf)}))
  
  WriteGenotype(AFRg, paste0(resdir_data, "AFRg.txt"))
  WriteGenotype(EURg, paste0(resdir_data, "EURg.txt"))
  
  X = cbind(1, Admprop, AFRg, EURg)
  
  res = rep(NA, 4)
  while (anyNA(res)){
    beta02use = order(sapply(beta0space, function(beta0){abs(mean(1/(1 + exp(- X %*% c(beta0, betas)))) - 0.2)}), decreasing = F)[1]
    allbetas = c(beta0space[beta02use], betas)
    
    
    probs = 1/(1 + exp(-mvrnorm(1, X %*% allbetas, tau2 * Kin)))
    pheno = sapply(probs, function(prob){rbinom(n = 1, size = 1, prob = prob)} )
    
    PHENO = data.frame(Admprop, pheno, ID = paste0("Sample", 1:n))
    
    model0_gmmat_true = glmmkin(fixed = pheno ~  1 + Admprop, data = PHENO, id = "ID", kins = Kin, family = binomial(link = "logit"))
    TractorMix.score(model0_gmmat_true, infilesTractorMix, outTractorMix)
    TractorMix_score_res = as.numeric(read.csv(outTractorMix, sep = "\t")[,10:11])
    
    Tractormix.wald_restbl = TractorMix.wald(infiles = infilesTractorMix, 
                                     data = PHENO, 
                                     fixed = pheno~ Admprop, kin = Kin, id = ID, family = binomial(link = "logit"))
    
    TractorMix_wald_res = c(Tractormix.wald_restbl$AFRg[,"P"], 
                            Tractormix.wald_restbl$EURg[,"P"])
      
    res = c(TractorMix_score_res, TractorMix_wald_res)
    names(res) = c("AFR_score", "EUR_score", "AFR_wald", "EUR_wald")
  }
  file.remove(paste0(resdir_data, "AFRg.txt"), paste0(resdir_data, "EURg.txt"), paste0(resdir_data, "TractorMix.res"))
  return(res)
}

pval_res_null_dich =  t(sapply(1:1000, RunSimulation_dich, betas_null_dich))
pval_res_alt_dich =  t(sapply(1:1000, RunSimulation_dich, betas_alt_dich))


########## writing results #############
write.table(pval_res_null_cont, paste0(resdir, "null_cont.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(pval_res_alt_cont, paste0(resdir, "alt_cont.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(pval_res_null_dich, paste0(resdir, "null_dich.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(pval_res_alt_dich, paste0(resdir, "alt_dich.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
