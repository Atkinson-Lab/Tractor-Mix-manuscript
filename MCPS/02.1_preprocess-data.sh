
#!/bin/bash
#SBATCH --job-name=preprocess-1-tracmix-navlg
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk

# Directories
base_dir=/well/emberson/users/pom143/projects/LAI_TractorMix
data_dir=$base_dir/geno-files/merged-hgdp-1kg-mcps
work_dir=/well/emberson/users/pom143/projects/LAI_TractorMix
out_dir=$work_dir/output-files
bcftools=/apps/well/bcftools/1.4.1/bin/bcftools
plink2=/well/emberson/shared/software/plink2/plink2

# Create reference panel sample list (HGDP, 1K Genome and MCPS)
awk '{print $2}' $data_dir/merged_hgdp-1kg-mcps_autosomes.fam > $out_dir/ref_mcps.txt

# Create reference panel sample list (HGDP and 1K Genome)
awk '!/^MCPS/' $out_dir/ref_mcps.txt > $out_dir/ref_1kg_hgdp.txt

# Create mcps panel sample list
awk '/^MCPS/' $out_dir/ref_mcps.txt > $out_dir/mcps_only.txt

# Create reference panel sample list (HGDP and 1K Genome)
awk 'NR==FNR{if($3 == "AMERICA" || $3 == "EUROPE") include[$1]; next} ($1 in include) {print $1}' $out_dir/reference-population-labels.txt $out_dir/ref_1kg_hgdp.txt > $out_dir/ref_subset_eur_amr.txt

# Reference to use in RFMIX and Tractor Mix 
awk 'NR==FNR{if($3 == "AMERICA" || $3 == "EUROPE") include[$1]=$3; next} ($1 in include) {print $1"\t"include[$1]}' $out_dir/reference-population-labels.txt $out_dir/ref_1kg_hgdp.txt > $out_dir/ref_label_eur_amr.txt

# Generate BMI sample list
awk '{print $1}' $base_dir/input-files/extracted-phenotypes-all-2.5perc-bmi.txt > $out_dir/ancestry_mcps_amr_eur_bmi.txt
cat $out_dir/mcps_ref_subset_eur_amr.txt $out_dir/ancestry_mcps_amr_eur_bmi.txt > $out_dir/refs_mcps_ancestry_amr_eur_bmi.txt

awk 'BEGIN {FS = OFS = "\t"} NR==1 {print "ID", "BMI", "AGE", "SEX", "Score_AMR"; next} {print $2, $3, $4, $5, $25}' $base_dir/input-files/extracted-phenotypes-all-2.5perc-bmi.txt > $out_dir/bmi.txt

