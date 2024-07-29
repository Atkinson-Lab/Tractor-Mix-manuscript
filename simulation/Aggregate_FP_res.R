args <- commandArgs(trailingOnly = TRUE)
inputtype <- args[1]
outputdir <- args[2]

library(data.table)
resfiles = list.files(path = paste0("/storage/atkinson/home/u241265/Tractor-Mix/simu/Sparse/", inputtype, "/"),recursive = TRUE, full.names =T)

Aggres <- function(resfile){
  res = fread(resfile, sep = "\t")
  res = res[complete.cases(res),]
  return(c(ntests = nrow(res), 
           p1 = sum(res$P < 0.05),
           p2 = sum(res$P < 0.0005),
           p3 = sum(res$P < 0.000005),
           p4 = sum(res$P < 0.00000005),
           p1_anc1 = sum(res$Pval_anc1 < 0.05),
           p2_anc1 = sum(res$Pval_anc1 < 0.0005),
           p3_anc1 = sum(res$Pval_anc1 < 0.000005),
           p4_anc1 = sum(res$Pval_anc1 < 0.00000005),
           p1_anc2 = sum(res$Pval_anc2 < 0.05),
           p2_anc2 = sum(res$Pval_anc2 < 0.0005),
           p3_anc2 = sum(res$Pval_anc2 < 0.000005),
           p4_anc2 = sum(res$Pval_anc2 < 0.00000005)
           ))
}

finalres = colSums(t(sapply(resfiles, Aggres)))
write.table(t(finalres), outputdir, quote = F, sep = "\t", row.names = F, col.names=T)
