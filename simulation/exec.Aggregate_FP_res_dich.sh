#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nodelist=mhgcp-c00,mhgcp-c01,mhgcp-c02,mhgcp-m00
#SBATCH --partition=gpu
#SBATCH --time=05:00:00
#SBATCH --job-name=Aggregare_FP_res_dich
#SBATCH --output=Aggregare_FP_res_dich.%j.out
#SBATCH --mem=32gb

echo "start"
/opt/R/4.1.1/bin/Rscript Aggregate_FP_res.R Dich Dich_FP_res.tsv
echo "finish"
