#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --job-name=Comp_ancP
#SBATCH --output=Comp_ancP.%j.out
#SBATCH --mem=32gb

echo "start"
/opt/R/4.1.1/bin/Rscript Comp_ancP.R
echo "finish"
