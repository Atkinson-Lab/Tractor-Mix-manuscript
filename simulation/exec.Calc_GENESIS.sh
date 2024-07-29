#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=mhgcp-c02
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --job-name=Calc_GENESIS
#SBATCH --output=Calc_GENESIS.%j.out
#SBATCH --mem=32gb

echo "start"
/opt/R/4.1.1/bin/Rscript Calc_GENESIS.R
echo "finish"
