#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=mhgcp-c00,mhgcp-c01,mhgcp-c02,mhgcp-m00
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
#SBATCH --job-name=UKBB_processing
#SBATCH --output=UKBB_processing.%j.out
#SBATCH --mem=32gb

mkdir ~/Tractor-Mix/UKBB/data/ADMIXTURE
mkdir ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB

echo "begain plink2 prune & remove duplicates ..."
/storage/atkinson/shared_software/software/plink2/plink2 \
    --vcf /storage/atkinson/shared_resources/datasets/UKB/QCedGenotypes/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto.vcf.bgz \
    --rm-dup exclude-all  \
    --set-all-var-ids @:#:\$r:\$a \
    --indep-pairwise 500kb 0.01 \
    --out ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto 


# extract common variants 

echo "intersect common variants from UKBB and TGP ... (a more careful analysis should take care of REF/ALT swap)"


