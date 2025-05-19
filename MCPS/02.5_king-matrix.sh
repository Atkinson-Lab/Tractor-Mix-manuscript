#!/bin/bash
#SBATCH --job-name=king
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=20:00:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk


# Modules
module purge                                             
module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate local_ancestry


# Define directories and input files
base_dir=/well/emberson/users/pom143/projects/LAI_TractorMix
input_dir=$base_dir/input-files
out_dir=$base_dir/output-files
king_out_dir=$out_dir/king

mkdir -p $king_out_dir

plink2=/well/emberson/shared/software/plink2/plink2
geno_dir=/well/emberson/projects/mcps/data/shortcuts/gsav2-chip-qcd
geno_file=$geno_dir/MCPS_Freeze_150.GT_hg38.pVCF.revised_qc.autosomes


# Generate Plink file set for analyses set 1 (n=10K, BMI)
#$plink2 --bfile $geno_file --keep $input_dir/extracted-phenotypes-all-2.5perc-bmi_renamed.txt --make-bed --out $out_dir/analysis-set-1_bmi

# King analysis for analysis set 1 
king -b $out_dir/analysis-set-1_bmi.bed --kinship --degree 3 --prefix $king_out_dir/analysis-set-1_bmi


#sed 's/FID.x/FID/g' extracted-phenotypes-all-2.5perc-bmi.txt > extracted-phenotypes-all-2.5perc-bmi_renamed.txt
#awk 'NR==1 {print "FID\tIID"} NR>1 {print $1 "\t" $1}' ancestry_mcps_amr_eur_bmi.txt > ancestry_mcps_amr_eur_bmi_with_IID.txt

