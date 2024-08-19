#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=mhgcp-c00,mhgcp-c01,mhgcp-c02,mhgcp-m00
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --job-name=TCP_processing
#SBATCH --output=TCP_processing.%j.out
#SBATCH --mem=32gb



echo "extract continental AFR and EUR from TGP"
mkdir ~/Tractor-Mix/UKBB/data/ADMIXTURE
mkdir ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP

unset PYTHONPATH
conda activate py37
# bcftools 1.20 installed

for CHR in {1..22}
do
   echo /storage/atkinson/shared_resources/reference/ReferencePanels/TGP/TGP_Phase3_GRCh37/FilteredVCFs/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.vcf.gz >> ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.filename
done


echo "concatenate chromosomes"
bcftools concat \
    --file-list ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.filename \
    --threads 16 \
    --output-type z\
    --output ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.vcf.gz

echo "reset var ID"
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.vcf.gz \
    --output-type z\
    --output ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.resetVarID.vcf.gz

echo "extract var ID"
bcftools query -f '%ID\n' ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.resetVarID.vcf.gz \
 > ~/Tractor-Mix/UKBB/data/ADMIXTURE/TGP/AFR_EUR_phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.resetVarID.vars