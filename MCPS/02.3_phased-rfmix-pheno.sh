#!/bin/bash
#SBATCH --job-name=phased-lai-navlg
#SBATCH --cpus-per-task=8
#SBATCH --time=18:00:00
#SBATCH --mem=128G
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk
#SBATCH --array=1-22

# Record the start time
start_time=$SECONDS

# Load necessary modules and activate the environment
module purge
module load Python/3.11.3-GCCcore-12.3.0
module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate local_ancestry

# Directories


gmap_dir_p=/well/emberson/users/bjk420/projects/popgen_v2/06_phasing/gmap_files
work_dir=/well/emberson/users/pom143/projects/LAI_TractorMix
inp_dir=$work_dir/input-files
out_dir=$work_dir/output-files
chr_out_auto=$out_dir/chr_autosomes_sample
chr_out_pha=$out_dir/chr_phased_sample
chr_out_rf=$out_dir/rfmix_lai_amr_eur
bcftools=/apps/well/bcftools/1.4.1/bin/bcftools
plink2=/well/emberson/shared/software/plink2/plink2

chr=${SLURM_ARRAY_TASK_ID}


process_data() {
    local trait=$1

    # Phase data
    SHAPEIT5_phase_common --input $chr_out_auto/modified_merged-refs-mcps_autosomes_${trait}_chr${chr}.vcf.gz \
             --map $gmap_dir_p/chr${chr}.b38.gmap.gz --region $chr --thread 16 \
             --output $chr_out_pha/ref_mcps_sample_amr_eur_phased_${trait}_chr${chr}.bcf \
             --log $chr_out_pha/logs/shapeit5_${trait}_chr${chr}.log

    # Compress and index the phased VCF (commented out here)
    $bcftools view $chr_out_pha/ref_mcps_sample_amr_eur_phased_${trait}_chr${chr}.bcf -O z -o $chr_out_pha/ref_mcps_sample_amr_eur_phased_${trait}_chr${chr}.vcf.gz
    $bcftools index $chr_out_pha/ref_mcps_sample_amr_eur_phased_${trait}_chr${chr}.vcf.gz &
    wait $!

    # Extract phased samples (commented out here)
    $bcftools view --force-samples -S $out_dir/ref_subset_eur_amr.txt $chr_out_pha/ref_mcps_sample_amr_eur_phased_${trait}_chr${chr}.vcf.gz -O z -o $chr_out_pha/ref_sample_amr_eur_phased_${trait}_chr${chr}.vcf.gz
    $bcftools view --force-samples -S $out_dir/mcps_${trait}_only.txt $chr_out_pha/ref_mcps_sample_amr_eur_phased_${trait}_chr${chr}.vcf.gz -O z -o $chr_out_pha/mcps_sample_amr_eur_phased_${trait}_chr${chr}.vcf.gz &
    wait $!

    # Check if RFMix succeeded
    if [ $? -eq 0 ]; then
        echo "RFMix completed successfully for chromosome $chr and trait $trait."

        # Run the first script and log the job ID
        job_1=$(sbatch --parsable -p short $work_dir/src/development/rfmix-amr-eur-bmi_vf.sh $chr $trait)

        # Check if the job was submitted successfully
        if [ $? -eq 0 ]; then
            echo "Submitted rfmix-amr-eur-bmi_vf.sh as job $job_1."
        else
            echo "Failed to submit rfmix-amr-eur-bmi_vf.sh for chromosome $chr and trait $trait."
        fi
    else
        echo "RFMix failed for chromosome $chr and trait $trait. Skipping dependent jobs."
    fi
}

for trait in "bmi"; do
    process_data $trait
done

# Calculate and print the duration for the current job
duration=$(( SECONDS - start_time ))
echo "The job for chromosome $chr took $((duration / 3600)) hours, $(((duration % 3600) / 60)) minutes, and $((duration % 60)) seconds."

