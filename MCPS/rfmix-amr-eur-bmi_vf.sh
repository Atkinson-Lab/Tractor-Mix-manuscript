#!/bin/bash
#SBATCH --job-name=rfmix-amr-eur-navlg
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk

# Modules
module purge
module load Python/3.11.3-GCCcore-12.3.0
module load Anaconda3/2024.02-1 
eval "$(conda shell.bash hook)"
conda activate local_ancestry

# Directories
gmap_dir=/well/emberson/users/bjk420/projects/popgen/04_rfmix/version_1/gmap_files
work_dir=/well/emberson/users/pom143/projects/LAI_TractorMix
inp_dir=$work_dir/input-files
out_dir=$work_dir/output-files
chr_out_auto=$out_dir/chr_autosomes_sample
chr_out_pha=$out_dir/chr_phased_sample
chr_out_rf=$out_dir/rfmix_lai_amr_eur


chr="$1"
trait="$2"

# Run RFMix for chromosomes 1 to 22

echo "Running RFMix for chromosome $chr and trait $trait..."
rfmix -f "$chr_out_pha/mcps_sample_amr_eur_phased_${trait}_chr${chr}.vcf.gz" \
      -r "$chr_out_pha/ref_sample_amr_eur_phased_${trait}_chr${chr}.vcf.gz" \
      --chromosome="$chr" \
      -m "$out_dir/ref_label_eur_amr.txt" \
      -g "$gmap_dir/chr${chr}.b38.reformat.gmap" \
      -n 2 \
      -G 15 \
      -num-threads 8 \
      --reanalyze-reference \
      -o "$chr_out_rf/local_admix_rfmix_amr_eur_${trait}_chr${chr}"


echo "All processing completed part rfmix for CHR_NUMBER $chr"
