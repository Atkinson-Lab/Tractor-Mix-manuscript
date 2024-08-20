#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=mhgcp
#SBATCH --nodelist=mhgcp-c00,mhgcp-c01,mhgcp-c02,mhgcp-m00
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --job-name=GMMAT
#SBATCH --output=GMMAT.%j.out
#SBATCH --mem=32gb
#SBATCH --array=1-22

echo "start GMMAT..."
/opt/R/4.1.1/bin/Rscript TC_LDL_SCA_GMMAT.R $SLURM_ARRAY_TASK_ID
echo "finish GMMAT..."