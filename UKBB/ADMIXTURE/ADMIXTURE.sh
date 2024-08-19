#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=mhgcp-c00,mhgcp-c01,mhgcp-c02,mhgcp-m00
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --time=06:30:00
#SBATCH --job-name=Run_ADMIXTURE
#SBATCH --output=Run_ADMIXTURE.%j.out
#SBATCH --mem=32gb


echo "PLINK extract common variants for UKBB"

/storage/atkinson/shared_software/software/plink2/plink2 \
    --vcf /storage/atkinson/shared_resources/datasets/UKB/QCedGenotypes/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto.vcf.bgz \
    --set-all-var-ids @:#:\$r:\$a \
    --extract ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB_TGP_intersect.vars \
    --export vcf bgz \
    --out ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-ldprune-intersect


echo "PLINK extract common variants for TGP"

/storage/atkinson/shared_software/software/plink2/plink2 \
    --vcf ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.vcf.gz \
    --extract ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB_TGP_intersect.vars \
    --export vcf bgz \
    --out ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.intersect


echo "sort and index UKBB data"
/storage/atkinson/shared_software/software/bcftools-1.19/bcftools \
    sort \
    ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-ldprune-intersect.vcf.gz \
    -Oz -o \
    ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-ldprune-intersect-sorted.vcf.gz

/storage/atkinson/shared_software/software/bcftools-1.19/bcftools \
    index -t \
    ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-ldprune-intersect-sorted.vcf.gz


echo "sort and index TGP data"

/storage/atkinson/shared_software/software/bcftools-1.19/bcftools \
    sort \
    ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.intersect.vcf.gz \
    -Oz -o \
    ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.intersect.sorted.vcf.gz

/storage/atkinson/shared_software/software/bcftools-1.19/bcftools \
    index -t \
    ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.intersect.sorted.vcf.gz



echo "merge TGP and UKBB"
/storage/atkinson/shared_software/software/bcftools-1.19/bcftools \
    merge ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-ldprune-intersect-sorted.vcf.gz \
    ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.intersect.sorted.vcf.gz \
    -o ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB_TGP_merged.vcf.gz


echo "convert merged vcf to plink1 format"
/storage/atkinson/shared_software/software/plink2/plink2 \
    --vcf ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB_TGP_merged.vcf.gz \
    --make-bed \
    --out ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB_TGP_merged


echo "run ADMIXTURE"

/storage/atkinson/shared_software/software/admixture_linux-1.3.0/admixture \
    ~/Tractor-Mix/UKBB/data/ADMIXTURE/UKBB_TGP_merged.bed \
    2 


