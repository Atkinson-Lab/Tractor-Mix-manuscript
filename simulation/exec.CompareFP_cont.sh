#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --job-name=CompareFP_cont
#SBATCH --output=CompareFP_cont.%j.out
#SBATCH --mem=32gb
#SBATCH --array=3-20

echo "start"
/opt/R/4.1.1/bin/Rscript CompareFP_cont.R $SLURM_ARRAY_TASK_ID 1
/opt/R/4.1.1/bin/Rscript CompareFP_cont.R $SLURM_ARRAY_TASK_ID 2
/opt/R/4.1.1/bin/Rscript CompareFP_cont.R $SLURM_ARRAY_TASK_ID 3
/opt/R/4.1.1/bin/Rscript CompareFP_cont.R $SLURM_ARRAY_TASK_ID 4
/opt/R/4.1.1/bin/Rscript CompareFP_cont.R $SLURM_ARRAY_TASK_ID 5
/opt/R/4.1.1/bin/Rscript CompareFP_cont.R $SLURM_ARRAY_TASK_ID 6
echo "finish"
