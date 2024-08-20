#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=mhgcp-c00,mhgcp-c01,mhgcp-c02,mhgcp-m00
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
#SBATCH --job-name=Make_King
#SBATCH --output=Make_King.%j.out
#SBATCH --mem=32gb

# compute King-robust, which will be later used for GENESIS 
# requires king, king.id 

mkdir ~/Tractor-Mix/UKBB/data/TRACTOR
mkdir ~/Tractor-Mix/UKBB/data/TRACTOR/GENESIS1

# 869,325 variants
# create a King matrix
# create a Plink bed files

/storage/atkinson/shared_software/software/plink2/plink2 \
        --vcf ~/Tractor-Mix/UKBB/data/PC_KING/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-prune.vcf.gz \
        --make-king square \
        --out ~/Tractor-Mix/UKBB/data/TRACTOR/GENESIS1/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-prune



/storage/atkinson/shared_software/software/plink2/plink2 \
        --vcf ~/Tractor-Mix/UKBB/data/PC_KING/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-prune.vcf.gz \
        --make-bed \
        --out ~/Tractor-Mix/UKBB/data/TRACTOR/GENESIS1/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-prune












echo "begain plink2 prune & remove duplicates ..."
/storage/atkinson/shared_software/software/plink2/plink2 \
    --vcf /storage/atkinson/shared_resources/datasets/UKB/QCedGenotypes/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto.vcf.bgz \
    --rm-dup exclude-all  \
    --indep-pairwise 500kb 0.2 \
    --out ~/Tractor-Mix/UKBB/data/TRACTOR/GENESIS1/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto 

echo "begain extract independent markers..."
/storage/atkinson/shared_software/software/plink2/plink2 \
    --vcf /storage/atkinson/shared_resources/datasets/UKB/QCedGenotypes/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto.vcf.bgz \
    --extract ~/Tractor-Mix/UKBB/data/TRACTOR/GENESIS1/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto.prune.in \
    --keep-allele-order \
    --export vcf bgz \
    --out ~/Tractor-Mix/UKBB/data/PC_KING/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-prune


echo "begain plink2 compute king robust..."
/storage/atkinson/shared_software/software/plink2/plink2 \
        --vcf ~/Tractor-Mix/UKBB/data/PC_KING/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-prune.vcf.gz \
        --make-king-table \
        --out ~/Tractor-Mix/UKBB/res/PC_KING/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-prune


