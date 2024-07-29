#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nodelist=mhgcp-g00
#SBATCH --partition=gpu
#SBATCH --job-name=Aggregare_FP_res_cont
#SBATCH --output=Aggregare_FP_res_cont.%j.out
#SBATCH --mem=32gb

echo "start"
/opt/R/4.1.1/bin/Rscript Aggregate_FP_res.R Cont Cont_FP_res.tsv
echo "finish"
