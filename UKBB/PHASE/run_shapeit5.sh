#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=mhgcp-c02
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --job-name=Shapeit5
#SBATCH --output=Shapeit5.%j.out
#SBATCH --mem=32gb
#SBATCH --array=1-22


echo "Set enviroment, initiate conda"
unset PYTHONPATH
~/miniconda3/bin/conda init bash
conda activate py37

export BCFTOOLS_PLUGINS=/storage/atkinson/shared_software/software/bcftools-1.19/plugins

echo "subset samples..."
/storage/atkinson/shared_software/software/bcftools-1.19/bcftools view -S ~/Tractor-Mix/UKBB/data/PHASE/TGP/TGP_AFR_EUR_4bcftools.tsv \
    /storage/atkinson/shared_resources/reference/ReferencePanels/TGP/TGP_Phase3_GRCh37/FilteredVCFs/ALL.chr$SLURM_ARRAY_TASK_ID.phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.vcf.gz \
    > ~/Tractor-Mix/UKBB/data/PHASE/TGP/AFR_EUR.chr$SLURM_ARRAY_TASK_ID.phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.vcf.gz


echo "compute AC for TGP..."
/storage/atkinson/shared_software/software/bcftools-1.19/bcftools +fill-tags ~/Tractor-Mix/UKBB/data/PHASE/TGP/AFR_EUR.chr$SLURM_ARRAY_TASK_ID.phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.vcf.gz \
    -Ob -o ~/Tractor-Mix/UKBB/data/PHASE/TGP/AFR_EUR.chr$SLURM_ARRAY_TASK_ID.phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.AC.vcf.gz -- -t AN,AC


echo "index TGP ..."
/storage/atkinson/shared_software/software/bcftools-1.19/bcftools \
    index -t \
    ~/Tractor-Mix/UKBB/data/PHASE/TGP/AFR_EUR.chr$SLURM_ARRAY_TASK_ID.phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.AC.vcf.gz



echo "compute AC for UKBB..."
/storage/atkinson/shared_software/software/bcftools-1.19/bcftools +fill-tags ~/Tractor-Mix/UKBB/data/PHASE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID.vcf.gz \
    -Ob -o ~/Tractor-Mix/UKBB/data/PHASE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID-AC.vcf.gz -- -t AN,AC


echo "index UKBB ..."
/storage/atkinson/shared_software/software/bcftools-1.19/bcftools \
    index -t \
    ~/Tractor-Mix/UKBB/data/PHASE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID-AC.vcf.gz


echo "run Shapeit5..."
SHAPEIT5_phase_common \
    --input ~/Tractor-Mix/UKBB/data/PHASE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID-AC.vcf.gz \
    --reference ~/Tractor-Mix/UKBB/data/PHASE/TGP/AFR_EUR.chr$SLURM_ARRAY_TASK_ID.phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.AC.vcf.gz\
    --region $SLURM_ARRAY_TASK_ID\
    --map ~/Tractor-Mix/UKBB/data/PHASE/UKBB/gmaps/chr$SLURM_ARRAY_TASK_ID.b37.gmap.gz \
    --output ~/Tractor-Mix/UKBB/data/PHASE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID-AC-phased.bcf \
    --thread 16

echo "convert to vcf.gz from bcf"
/storage/atkinson/shared_software/software/bcftools-1.19/bcftools \
    convert -O z -o ~/Tractor-Mix/UKBB/data/PHASE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID-AC-phased.vcf.gz \
    ~/Tractor-Mix/UKBB/data/PHASE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR$SLURM_ARRAY_TASK_ID-AC-phased.bcf


echo "finish..."
