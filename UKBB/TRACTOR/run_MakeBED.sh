#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=mhgcp-c02
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --time=1:30:00
#SBATCH --job-name=Make_Bed
#SBATCH --output=Make_Bed.%j.out
#SBATCH --mem=16gb
#SBATCH --array=1-22


unset PYTHONPATH
~/miniconda3/bin/conda init bash
conda activate py37


/storage/atkinson/shared_software/software/plink2/plink2 \
    --vcf ~/Tractor-Mix/UKBB/data/PHASE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID-AC-phased.vcf.gz \
    --make-bed \
    --out ~/Tractor-Mix/UKBB/data/TRACTOR/BEDs/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID-AC-phased

