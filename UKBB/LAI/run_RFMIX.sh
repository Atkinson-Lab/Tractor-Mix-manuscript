#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=mhgcp-c02
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --job-name=RFMIX2
#SBATCH --output=RFMIX2.%j.out
#SBATCH --mem=32gb
#SBATCH --array=1-22

echo "Run rfmix2..."
/storage/atkinson/shared_software/software/rfmix/rfmix \
    -f ~/Tractor-Mix/UKBB/data/PHASE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID-AC-phased.vcf.gz \
    -r ~/Tractor-Mix/UKBB/data/PHASE/TGP/AFR_EUR.chr$SLURM_ARRAY_TASK_ID.phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.AC.vcf.gz \
    -m ~/Tractor-Mix/UKBB/data/LAI/TGP_AFR_EUR.tsv \
    -g ~/Tractor-Mix/UKBB/data/LAI/gmaps_rfmix/chr$SLURM_ARRAY_TASK_ID.b37.gmap \
    -o ~/Tractor-Mix/UKBB/data/LAI/msp/UKBB-AFR-deconvoluted-chr$SLURM_ARRAY_TASK_ID \
    --chromosome=$SLURM_ARRAY_TASK_ID

echo "finish..."
