#!/bin/bash
#SBATCH --job-name=dosages-navlg
#SBATCH --time=00:30:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk
#SBATCH --array=1-7

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
    
    echo "Submitted dosages."

    # Run the second script after the first script has finished
    sbatch -p short $work_dir/src/development/02.4_extract-tracts.sh $chr $trait

}

for trait in "bmi"; do
    process_data $trait
done

# Calculate and print the duration for the current job
duration=$(( SECONDS - start_time ))
echo "The job for chromosome $chr took $((duration / 3600)) hours, $(((duration % 3600) / 60)) minutes, and $((duration % 60)) seconds."

