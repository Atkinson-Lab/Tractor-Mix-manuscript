#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --job-name=Extract_Tracts
#SBATCH --output=Extract_Tracts.%j.out
#SBATCH --mem=32gb
#SBATCH --array=1-22


unset PYTHONPATH
~/miniconda3/bin/conda init bash
conda activate py37

echo "extract tracts..."
python ~/Tractor-Mix/UKBB/src/Tractor/scripts/ExtractTracts.py \
      --msp ~/Tractor-Mix/UKBB/data/LAI/msp/UKBB-AFR-deconvoluted-chr$SLURM_ARRAY_TASK_ID \
      --vcf-prefix ~/Tractor-Mix/UKBB/data/PHASE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID-AC-phased\
      --zipped \
      --num-ancs 2 \
      --output-path  ~/Tractor-Mix/UKBB/data/TRACTOR/HAP_DOSE/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID-AC-phased

echo "finish..."
