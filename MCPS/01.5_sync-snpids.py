#!/usr/bin/python -O
# Jason Matthew Torres & Alejandra Vergara-Lope
'''

module load Python/3.11.3-GCCcore-12.3.0
python 01.4_sync-snpids.py

'''
# libraries
import sys,os
from shutil import copyfile
import subprocess as sp

work_dir = "/well/emberson/users/pom143/projects/LAI_TractorMix/"
plink2 = "/well/emberson/shared/software/plink2/plink2"

def sync_snps(study_dir,study_name):
    sys.stdout.write("\nRunning dataset: %s\n" % study_name)
    bim_file = study_dir + study_name + ".bim"
    copyfile(bim_file,bim_file+"_copy")
    fin = open(bim_file+"_copy",'r')
    fout = open(bim_file+"_temp",'w')
    fout2 = open(study_dir + study_name + ".indels.txt",'w')
    count = 0
    for line in fin:
        count += 1
        sys.stdout.write("\r")
        sys.stdout.write("Count: %d" % count)
        sys.stdout.flush()
        l = line.strip().split()
        chrom,pos = l[0],l[3]
        a1,a2 = l[4],l[5]
        a_list = [len(a1),len(a2)]
        a_list.sort(reverse=True)
        max_a_num = a_list[0]
        snpid = chrom+":"+pos
        if max_a_num > 1:
            snpid = snpid + ":"+str(max_a_num)
            fout2.write(snpid+"\n")
        l[1] = snpid
        fout.write("\t".join(l)+"\n")
    fin.close()
    fout.close()
    fout2.close()
    os.rename(bim_file+"_temp",bim_file)

def plink_subset(study_dir,study_name):
    prefix = study_dir+study_name
    indel_file = study_dir+study_name+".indels.txt"
    command = [plink2, "--bfile", prefix, "--exclude", indel_file, "--max-alleles", "2",\
    "--make-bed","--out",prefix+".biallelic"]
    sp.check_call(command)


def main():
    study_dir1 = work_dir + "geno-files/merged-hgdp-1kg/"
    study_dir2 = work_dir + "geno-files/mcps/"
    ## Autosomes
    sys.stdout.write("\nAutosome variants...\n")
    study_name1 = "merged_hgdp-1kg_autosomes"
    study_name2 = "topmed-imputed-genotypes-2.5P-ND"
    sync_snps(study_dir1,study_name1)
    sync_snps(study_dir2,study_name2)
    plink_subset(study_dir1,study_name1)
    plink_subset(study_dir2,study_name2)

if (__name__=="__main__"):
     main()
