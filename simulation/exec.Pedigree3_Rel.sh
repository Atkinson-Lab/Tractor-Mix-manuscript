#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=mhgcp-c02
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --job-name=Admixture_Ind
#SBATCH --output=Admixture_Ind.%j.out
#SBATCH --mem=32gb
#SBATCH --array=26-50

echo "start"
/opt/R/4.1.1/bin/Rscript Pedigree3_Rel.R $SLURM_ARRAY_TASK_ID
echo "finish"
