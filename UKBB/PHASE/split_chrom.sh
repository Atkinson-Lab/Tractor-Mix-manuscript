#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --nodelist=mhgcp-c02
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --job-name=Split_vcf
#SBATCH --output=Split_vcf.%j.out
#SBATCH --mem=32gb

# echo "index the UKBB file"
# /storage/atkinson/shared_software/software/bcftools-1.19/bcftools index \
#     -t /storage/atkinson/shared_resources/datasets/UKB/QCedGenotypes/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto.vcf.bgz

echo "split UKBB into 22 vcfs"
/storage/atkinson/shared_software/software/bcftools-1.19/bcftools \
    index -s /storage/atkinson/shared_resources/datasets/UKB/QCedGenotypes/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto.vcf.bgz \
    | cut -f 1 | \
    while read C; do bcftools view -O z -o ~/Tractor-Mix/UKBB/data/PHASE/UKBB/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-CHR${C}.vcf.gz /storage/atkinson/shared_resources/datasets/UKB/QCedGenotypes/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto.vcf.bgz "${C}" ; done
