#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=mhgcp
#SBATCH --nodelist=mhgcp-c00,mhgcp-c01,mhgcp-c02,mhgcp-m00
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --job-name=make_PC_GRM_from_GENESIS
#SBATCH --output=make_PC_GRM_from_GENESIS.%j.out
#SBATCH --mem=32gb

echo "start..."
/opt/R/4.1.1/bin/Rscript ~/Tractor-Mix/UKBB/src/TRACTOR/make_PC_GRM_from_GENESIS.R
echo "finish..."
