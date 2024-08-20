library(data.table)
library(dplyr)
hapdose_id = read.csv("~/Tractor-Mix/UKBB/data/TRACTOR/PHE/hapdose_id.txt", header = F)
colnames(hapdose_id) = c("sample_id")
biomarker = fread("~/Tractor-Mix/UKBB/data/TRACTOR/PHE/Phenotypes_Everyone_uk_round2_allSamples_biomarkers_phesant_QC.tsv", sep = "\t")
phecode_admix = read.csv("~/Tractor-Mix/UKBB/data/TRACTOR/PHE/phecode_admix_wide.tsv", sep = "\t", check.names = FALSE)
phe1 = fread("~/Tractor-Mix/UKBB/data/TRACTOR/PHE/Phenotypes_Everyone_uk_round2_allSamples_phenos_phesant_QC.1.tsv", sep = "\t")

Phenotypes2use = left_join(hapdose_id, biomarker, by = c("sample_id" = "userId")) %>% 
                    left_join(phecode_admix, by = c("sample_id" = "userId")) %>% 
                    left_join(phe1, by = c("sample_id" = "userId")) %>%
                    select(c("sample_id", "age.y", "sex.y","30690", "30780", "282.5")) %>% 
                    rename("age" = "age.y", "sex" = "sex.y", "TC" = "30690", "LDL" = "30780", "SCA" = "282.5") %>% 
                    mutate("sample_id" = as.character(sample_id), "SCA" = as.integer(SCA * 1), 
                           "TC" = as.numeric(TC), LDL = as.numeric(LDL))

load("~/Tractor-Mix/UKBB/data/TRACTOR/GENESIS/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-prune-PC_GRMs.RData")
PCAirmat = as.data.frame(PC_GRMs$PCAirmat)
PCAirmat$sample_id = row.names(PCAirmat)
PCRelatemat = PC_GRMs$PCRelatemat



# ensure same order of Phenotype, PC, and Kinship 
all.equal(Phenotypes2use$sample_id, PCAirmat$sample_id, 
          row.names(PCRelatemat), colnames(PCRelatemat))


### merge PCAir with 3 phenotypes 
TC_PCAir = inner_join(Phenotypes2use, PCAirmat, by = c("sample_id" = "sample_id"))  %>% 
  select("sample_id", "TC", "age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5")
LDL_PCAir = inner_join(Phenotypes2use, PCAirmat, by = c("sample_id" = "sample_id"))  %>% 
  select("sample_id", "LDL", "age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5")
SCA_PCAir = inner_join(Phenotypes2use, PCAirmat, by = c("sample_id" = "sample_id"))  %>% 
  select("sample_id", "SCA", "age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5")




library(GMMAT)

### Fit the null model for 3 phenotype, use PCAir + PCRelate
TC_Null = glmmkin(fixed = TC ~  age + sex + PC1 + PC2 + PC3 + PC4 + PC5, 
                           data = TC_PCAir, id = "sample_id", kins = PCRelatemat, 
                           family = gaussian(link = "identity"))
LDL_Null = glmmkin(fixed = LDL ~  age + sex + PC1 + PC2 + PC3 + PC4 + PC5, 
                           data = LDL_PCAir, id = "sample_id", kins = PCRelatemat, 
                           family = gaussian(link = "identity"))
SCA_Null = glmmkin(fixed = SCA ~  age + sex + PC1 + PC2 + PC3 + PC4 + PC5, 
                           data = SCA_PCAir, id = "sample_id", kins = PCRelatemat, 
                           family = binomial(link = "logit"))



save(TC_Null, file = "~/Tractor-Mix/UKBB/data/TRACTOR/NULLs/TC_Null.RData")
save(LDL_Null, file = "~/Tractor-Mix/UKBB/data/TRACTOR/NULLs/LDL_Null.RData")
save(SCA_Null, file = "~/Tractor-Mix/UKBB/data/TRACTOR/NULLs/SCA_Null.RData")

