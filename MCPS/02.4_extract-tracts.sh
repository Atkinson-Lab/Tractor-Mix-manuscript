#!/bin/bash
#SBATCH --job-name=ext1-2-navlg
#SBATCH -c 20
#SBATCH -p long
#SBATCH --mem=450G
#SBATCH --array=1-2

# Record the start time
start_time=$SECONDS

# Load necessary modules and activate the environment
module purge
module load Python/3.11.3-GCCcore-12.3.0
module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate local_ancestry

chr=${SLURM_ARRAY_TASK_ID}

base_dir=/well/emberson/users/pom143/projects/LAI_TractorMix
data_dir=$base_dir/geno-files/merged-hgdp-1kg-mcps
inp_dir=$base_dir/input-files
out_dir=$base_dir/output-files
chr_out_auto=$out_dir/chr_autosomes_sample
chr_out_pha=$out_dir/chr_phased_sample
chr_out_rf=$out_dir/rfmix_lai_amr_eur
vcffile=$chr_out_pha/mcps_sample_amr_eur_phased_bmi_chr${chr}.vcf.gz
mspfile=$chr_out_rf/local_admix_rfmix_amr_eur_bmi_chr${chr}.msp.tsv

python $base_dir/src/development/extract_tracts.py \
    --vcf $vcffile \
    --msp $mspfile \
    --num-ancs 2 \
    --output-dir ${chr_out_rf}
