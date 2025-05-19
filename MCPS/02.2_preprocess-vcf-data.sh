#!/bin/bash
#SBATCH --job-name=preprocess-2-data-vcf-tracmix-navlg
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=20:00:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk
#SBATCH --array=1-21

# Directories
work_dir=/well/emberson/users/pom143/projects/LAI_TractorMix
data_dir=$work_dir/geno-files/merged-hgdp-1kg-mcps
inp_dir=$work_dir/input-files
out_dir=$work_dir/output-files
chr_out_auto=$out_dir/chr_autosomes_sample
bcftools=/apps/well/bcftools/1.4.1/bin/bcftools
plink2=/well/emberson/shared/software/plink2/plink2


chr=${SLURM_ARRAY_TASK_ID}


$plink2 --bfile $data_dir/merged_hgdp-1kg-mcps_autosomes --chr ${chr} \
        --export vcf bgz id-delim="-" \
        --out $chr_out_auto/merged-refs-mcps_autosomes_chr${chr}

$bcftools index $chr_out_auto/merged-refs-mcps_autosomes_chr${chr}.vcf.gz

# Reheader the VCF file: example MCPS_MCPSRGN006190_MEXB006573-MCPS_MCPSRGN006190_MEXB006573 → MCPS_MCPSRGN006190_MEXB006573
$bcftools reheader -s <($bcftools query -l $chr_out_auto/merged-refs-mcps_autosomes_chr${chr}.vcf.gz | \
    sed -E -e 's/^0-//' \
    -e 's/(MCPS_[^_]+_[^_]+)-.*/\1/' | \
    awk '!seen[$0]++' | grep -v '^$') \
    $chr_out_auto/merged-refs-mcps_autosomes_chr${chr}.vcf.gz -o $chr_out_auto/modified_merged-refs-mcps_autosomes_chr${chr}.vcf

for trait in "bmi"; do

    $bcftools view --force-samples -S $out_dir/refs_mcps_ancestry_amr_eur_${trait}.txt \
    $chr_out_auto/modified_merged-refs-mcps_autosomes_chr${chr}.vcf \
     -O z -o $chr_out_auto/modified_merged-refs-mcps_autosomes_${trait}_chr${chr}.vcf.gz

    $bcftools index $chr_out_auto/modified_merged-refs-mcps_autosomes_${trait}_chr${chr}.vcf.gz
done

