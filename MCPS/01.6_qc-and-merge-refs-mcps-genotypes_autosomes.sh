#!/bin/bash
#SBATCH --job-name=merged-mcps-ref
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=20:00:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk

module load Python/3.11.3-GCCcore-12.3.0
plink=/well/emberson/shared/software/plink2/plink2
plink1=/well/emberson/shared/software/plink1/plink


### Create sub directories
workdir=/well/emberson/users/pom143/projects/LAI_TractorMix
basedir=/well/emberson/users/bjk420/projects/LAI_TractorMix
datadir1=$basedir/geno-files/merged-hgdp-1kg
datadir2=$workdir/geno-files/mcps

mkdir -p $workdir/geno-files/merged-hgdp-1kg-mcps
mergedir=$workdir/geno-files/merged-hgdp-1kg-mcps

name1='merged_hgdp-1kg_autosomes.biallelic'
name2='topmed-imputed-genotypes-2.5P-ND.biallelic'

### Run python script that syncs SNP ids between study bim files
### this file also now filter out indels and generates biallelic plink files
### note: this is only necessary because vcf file for build 38 1Kgenomes only has "." in the ID column
#python $workdir/01.2_sync-snpids.py


### Filter reference and study data for non A-T or G-C SNPs
echo 'Filter reference and study data for non AT or GC SNPs'
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $datadir1/$name1.bim > \
    $mergedir/$name1.ac_gt_snps

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $datadir2/$name2.bim > \
    $mergedir/$name2.ac_gt_snps

mkdir -p $mergedir/plink_log/

$plink --bfile $datadir1/$name1 \
      --exclude $mergedir/$name1.ac_gt_snps \
      --make-bed \
      --out $mergedir/$name1.no_ac_gt_snps
mv $mergedir/$name1.no_ac_gt_snps.log $mergedir/plink_log/$name1.no_ac_gt_snps.log

$plink --bfile $datadir2/$name2 \
      --exclude $mergedir/$name2.ac_gt_snps \
      --make-bed \
      --out $mergedir/$name2.no_ac_gt_snps
mv $mergedir/$name2.no_ac_gt_snps.log $mergedir/plink_log/$name2.no_ac_gt_snps.log


### Filter reference data for the same SNP set as in study
echo 'Filter reference data for the same SNP set as in study'
$plink --bfile $datadir2/$name2 \
      --extract $mergedir/$name1.no_ac_gt_snps.bim \
      --make-bed \
      --out $mergedir/$name2.pruned
mv $mergedir/$name2.pruned.log $mergedir/plink_log/$name2.pruned.log

$plink --bfile $datadir1/$name1 \
      --extract $mergedir/$name2.no_ac_gt_snps.bim \
      --make-bed \
      --out $mergedir/$name1.pruned
mv $mergedir/$name1.pruned.log $mergedir/plink_log/$name1.pruned.log


### Position mismatch
echo 'Position mismatch evaluate'
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4) {print a[$2],$2}' \
    $mergedir/$name1.pruned.bim $mergedir/$name2.pruned.bim > \
    $mergedir/${name2}.toUpdatePos

### Possible allele flips
echo 'Evaluate possible allele flips'
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $mergedir/$name1.pruned.bim $mergedir/$name2.pruned.bim > \
    $mergedir/$name2.toFlip

### There are duplicate variant IDs remaining in MCPS dataset
### Need to first remove these variants
echo 'Detecting and removing duplicates'
python $workdir/src/development/01.4_find-duplicate-variants.py $mergedir/$name2.pruned.bim $mergedir/$name2.duplicate-variants
$plink1 --bfile $mergedir/$name2.pruned \
       --exclude $mergedir/$name2.duplicate-variants \
       --make-bed \
       --out $mergedir/$name2.pruned.noDups
mv $mergedir/$name2.pruned.noDups.log $mergedir/plink_log/$name2.pruned.noDups.log


### Update positions and flip alleles
echo 'Update positions and flip alleles'
###--update-map $qcdir/$refname.toUpdatePos 1 2 \ # Will not run this arg
$plink1 --bfile $mergedir/$name2.pruned.noDups \
       --flip $mergedir/$name2.toFlip \
       --make-bed \
       --out $mergedir/$name2.flipped
mv $mergedir/$name2.flipped.log $mergedir/plink_log/$name2.flipped.log

### Remove mismatches
echo 'Remove mismatches'
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $mergedir/$name1.pruned.bim $mergedir/$name2.flipped.bim > \
    $mergedir/$name2.mismatch

$plink --bfile $mergedir/$name2.flipped \
       --exclude $mergedir/$name2.mismatch \
       --make-bed \
       --out $mergedir/$name2.clean
mv $mergedir/$name2.clean.log $mergedir/plink_log/$name2.clean.log

cat $mergedir/$name2.clean.bim | cut -f 2 | uniq -d > $mergedir/$name2.dup.snps

$plink --bfile $mergedir/$name2.clean \
       --exclude $mergedir/$name2.dup.snps \
       --make-bed \
       --out $mergedir/$name2.cleaner
mv $mergedir/$name2.cleaner.log $mergedir/plink_log/$name2.cleaner.log

### Merge study genotypes and reference data
echo 'Merge study genotypes and reference data'
cat $mergedir/$name2.cleaner.bim | cut -f 2 > $mergedir/$name2.cleaner.snps
$plink --bfile $mergedir/$name1.pruned \
       --extract $mergedir/$name2.cleaner.snps \
       --make-bed \
       --out $mergedir/$name1.cleaner

$plink1 --bfile $mergedir/$name1.cleaner \
      --bmerge $mergedir/$name2.cleaner.bed $mergedir/$name2.cleaner.bim \
        $mergedir/$name2.cleaner.fam \
      --make-bed \
      --out $mergedir/merged_hgdp-1kg-mcps_autosomes
#mv $mergedir/$name1.merge.$name2.log $mergedir/plink_log
