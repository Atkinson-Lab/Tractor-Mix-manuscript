#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --job-name=Split_convert_bed
#SBATCH --output=Split_convert_bed.%j.out
#SBATCH --mem=32gb

echo "convert to bed..."
/storage/atkinson/shared_software/software/plink2/plink2 \
    --vcf /storage/atkinson/shared_resources/datasets/UKB/QCedGenotypes/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto.vcf.bgz \
    --make-bed \
    --out ~/Tractor-Mix/UKBB/data/TRACTOR/BED/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto

echo "split chromosomes..."
for chr in {1..23}; do \
    /storage/atkinson/shared_software/software/plink2/plink2 \
    --bfile ~/Tractor-Mix/UKBB/data/TRACTOR/BED/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto \
    --chr $chr \
    --make-bed \
    --out ~/Tractor-Mix/UKBB/data/TRACTOR/BED/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto-chr${chr}; \
done

echo "finish"
