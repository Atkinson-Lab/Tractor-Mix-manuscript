# Scripts to generate samples by gene dropping, in the context of multi-ethenic and admixed population
# A person is defined as a list with 4 matricies: hap1, hap2, la1, la2
# hap1 and hap2 represents two haploids
# la1 and la2 represents2 two local ancestries

GenerateChild <- function(parent1, parent2){
  
  genolen = length(parent1$hap1)
  # recombination of parent1
  idx1 = sample(c(TRUE, FALSE), size = genolen, replace = T)
  childhap1 = rep(NA, genolen)
  childhap1[idx1] = parent1$hap1[idx1]
  childhap1[!idx1] = parent1$hap2[!idx1]
  childla1 = rep(NA, genolen)
  childla1[idx1] = parent1$la1[idx1]
  childla1[!idx1] = parent1$la2[!idx1]
  
  # recombination of parent2
  idx2 = sample(c(TRUE, FALSE), size = genolen, replace = T)
  childhap2 = rep(NA, genolen)
  childhap2[idx2] = parent2$hap1[idx2]
  childhap2[!idx2] = parent2$hap2[!idx2]
  childla2 = rep(NA, genolen)
  childla2[idx2] = parent2$la1[idx2]
  childla2[!idx2] = parent2$la2[!idx2]
  
  return(list(
    hap1 = childhap1,
    hap2 = childhap2,
    la1 = childla1,
    la2 = childla2
  ))
}


# generate minor allele frequency according to fst
GenerateMafs <- function(p, fst, runif1, runif2){
  f = runif(p, runif1, runif2)
  maf1 = rbeta(p, (1-fst)/fst*f, (1-fst)/fst*(1-f))
  maf2 = rbeta(p, (1-fst)/fst*f, (1-fst)/fst*(1-f))
  return(list(maf1 = maf1, maf2 = maf2))
}


## Generate independent Homogeneous population
GenerateHomo <- function(mafs, la){
  return(list(hap1 = sapply(mafs, function(maf){return(rbinom(1, 1, maf))}),
              hap2 = sapply(mafs, function(maf){return(rbinom(1, 1, maf))}),
              la1 = rep(la, length(mafs)),
              la2 = rep(la, length(mafs))))
}




## Generate independent Admixture population 
GenerateAdm <- function(mafs1, mafs2, prop){
  # generate local ancestry 
  # local ancestry = 0 is inherited from ancestry 1
  # local ancestry = 1 is inherited from ancestry 2
  la1 = rbinom(length(mafs1), size = 1, prob = (1 - prop))
  la2 = rbinom(length(mafs1), size = 1, prob = (1 - prop))
  
  # create empty haplotypes
  hap1 = rep(NA, length(mafs1))
  hap2= rep(NA, length(mafs1))
  
  # when local ancestry = 0, then use maf1 to generate genotype
  # when local ancestry = 1, then use maf2 to generate genotype
  hap1[which(la1 == 0)] = sapply(mafs1[which(la1 == 0)], function(prob){rbinom(n = 1, size = 1, prob)})
  hap1[which(la1 == 1)] = sapply(mafs2[which(la1 == 1)], function(prob){rbinom(n = 1, size = 1, prob)})
  
  hap2[which(la2 == 0)] = sapply(mafs1[which(la2 == 0)], function(prob){rbinom(n = 1, size = 1, prob)})
  hap2[which(la2 == 1)] = sapply(mafs2[which(la2 == 1)], function(prob){rbinom(n = 1, size = 1, prob)})
  
  return(list(
    hap1 = hap1,
    hap2 = hap2,
    la1 = la1,
    la2 = la2
  ))
}

# write genotype files in Tractor form
WriteGenotype <- function(G, filename){
  row.names(G) = paste0("Sample", 1:nrow(G))
  write.table(x = cbind(CHR = 1, POS  = 1:ncol(G), SNP = ".",  REF = "A", ALT = "G", t(G)),
              file = filename,
              quote = F, 
              sep = "\t",
              row.names = F,
              col.names = T)
}



# plot PCs 
# PlotPC <- function(G,Admprop){
#   PCs = prcomp(G)$x[,1:2]
#   df2plt = data.frame(cbind(PCs, AdmProp = Admprop))
#   p1 = ggplot(df2plt) + aes(x = PC1, y = PC2, colour = AdmProp) + geom_point() + theme_classic() + ggtitle("PC1 vs PC2")
#   p2 = ggplot(df2plt) + aes(x = PC1, y = AdmProp) + geom_point() + theme_classic() + ggtitle("PC1 vs Admixture")
#   return(p1 + p2)
# }


# lambda GC calculator
CalcLambdaGC = function(ps, df){
  expect_stats = qchisq(ppoints(length(ps)), df = df, lower = F)
  obs_stats = qchisq((ps), df = df, lower = F)
  return(list(expect_stats = expect_stats, obs_stats = obs_stats))
}


# seperate a genotype matrix into multiple pieces for testing and constructing GRM 
SeperateGenotype <- function(Genotype, Admprop, pGRM, filename){
  Geno4GRM = Genotype[, 1:pGRM]
  Geno4Test = Genotype[, (pGRM + 1):ncol(Genotype)]
  WriteGenotype(Geno4GRM,  paste0(filename, "_4GRM.tsv"))
  WriteGenotype(Geno4Test, paste0(filename, "_4Test.tsv"))
  return(list(Geno4GRM = Geno4GRM,
              Geno4Test = Geno4Test,
              Admprop = Admprop))
}




PrepareGDS <- function(Geno4GRM, filename){
  bim = as.data.frame(cbind(chr = rep(1, ncol(Geno4GRM)), id = ".", posg = 1:ncol(Geno4GRM), pos = 1:ncol(Geno4GRM), ref = "C", alt = "G"))
  fam = as.data.frame(cbind(fam = paste0("Sample", 1:nrow(Geno4GRM)), id = paste0("Sample", 1:nrow(Geno4GRM)), pat = ".", mat = ".", sex = ".", pheno = "." ))
  write_plink(filename, t(Geno4GRM), bim = bim, fam = fam)
  
  # create gds file
  snpgdsBED2GDS(bed.fn = paste0(filename, ".bed"), 
                bim.fn = paste0(filename, ".bim"), 
                fam.fn = paste0(filename, ".fam"), 
                out.gdsfn = paste0(filename, ".gds"))
  
  
  gds <- snpgdsOpen(paste0(filename, ".gds"))
  KINGmat = snpgdsIBDKING(gds)$kinship
  snpgdsClose(gds)
  
  
  gds <- GdsGenotypeReader(filename = paste0(filename, ".gds"))
  gdsData <- GenotypeData(gds)
  row.names(KINGmat) = paste0("Sample", 1:nrow(Geno4GRM))
  colnames(KINGmat) = paste0("Sample", 1:nrow(Geno4GRM))
  
  # PC Air calculation 
  PCair <- pcair(gdsData, kinobj = KINGmat, divobj = KINGmat)
  PCairpart <- pcairPartition(kinobj = KINGmat, divobj = KINGmat)
  
  # PC relate calculation 
  gdsData <- GenotypeBlockIterator(gdsData)
  PCrelate <- pcrelate(gdsData, pcs = PCair$vectors[,1:2], 
                          training.set = PCairpart$unrels,
                          BPPARAM = BiocParallel::SerialParam())
  
  
  GRM_PC = pcrelateToMatrix(PCrelate, sample.include = paste0("Sample", 1:nrow(Geno4GRM)), thresh = 2^(-11/2), scaleKin = 2)
  GRM_PC = round(as.matrix(GRM_PC)[paste0("Sample", 1:nrow(Geno4GRM)), paste0("Sample", 1:nrow(Geno4GRM))], 4)
  
  return(list(PCAir = PCair$vectors, 
              PCRelate = GRM_PC))
}



# get ancestry specific allele count
GetAncestrySpecCount <- function(ind, query){
    return(((ind$la1 == query) & (ind$hap1 == 1)) * 1 + ((ind$la2 == query) & (ind$hap2 == 1)) * 1)
}


# get global ancestry proportion with respect to population 0
GetGlobAncestry <- function(ind){
  return(mean(c(ind$la1, ind$la2) == 0))
}











