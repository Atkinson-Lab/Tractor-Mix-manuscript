library(data.table)

TGP_vars = "~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.resetVarID.vars"
UKBB_vars = "~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto.prune.in"

TGP_vars = fread(TGP_vars, header = FALSE) # 15,414,515
UKBB_vars = fread(UKBB_vars, header = FALSE) # 36,660
intersect_vars = intersect(TGP_vars$V1, UKBB_vars$V1) # 32,132

write.table(intersect_vars, "~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB_TGP_intersect.vars", col.names = F, row.names = F, quote = F)