library(data.table)

kin0 = fread("~/Tractor-Mix/UKBB/res/PC_KING/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-prune.kin0", sep = "\t")
# 1st degree relatedness 
sum(kin0$KINSHIP > 0.177)

# 2nd degree relatedness 
sum(kin0$KINSHIP > 0.088)