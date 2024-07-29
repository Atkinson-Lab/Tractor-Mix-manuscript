## This is a pedigree using PC-Air paper setup


MakePedigree <- function(MAFs, param1){
  Sample1 = GenerateAdm(MAFs$maf1, MAFs$maf2,  rbeta(1, param1, 20 - param1))
  Sample2 = GenerateAdm(MAFs$maf1, MAFs$maf2,  rbeta(1, param1, 20 - param1))
  Sample3 = GenerateAdm(MAFs$maf1, MAFs$maf2,  rbeta(1, param1, 20 - param1))
  Sample4 = GenerateChild(Sample1, Sample2)
  Sample5 = GenerateChild(Sample1, Sample2)
  Sample6 = GenerateAdm(MAFs$maf1, MAFs$maf2,  rbeta(1, param1, 20 - param1))
  Sample7 = GenerateChild(Sample1, Sample2)
  Sample8 = GenerateAdm(MAFs$maf1, MAFs$maf2,  rbeta(1, param1, 20 - param1))
  Sample9 = GenerateChild(Sample3, Sample4)
  Sample10= GenerateChild(Sample3, Sample4) 
  Sample11= GenerateChild(Sample3, Sample4) 
  Sample12= GenerateChild(Sample5, Sample6) 
  Sample13= GenerateChild(Sample5, Sample6) 
  Sample14= GenerateChild(Sample7, Sample8) 
  Sample15= GenerateChild(Sample7, Sample8) 
  Sample16= GenerateChild(Sample7, Sample8) 
  Sample17= GenerateAdm(MAFs$maf1, MAFs$maf2,  rbeta(1, param1, 20 - param1))
  Sample18= GenerateChild(Sample16, Sample17) 
  Sample19= GenerateChild(Sample16, Sample17) 
  Sample20= GenerateChild(Sample16, Sample17) 
  
  
  GenoMatTot = rbind(Sample9$hap1 + Sample9$hap2,
                      Sample10$hap1 + Sample10$hap2,
                      Sample11$hap1 + Sample11$hap2,
                      Sample12$hap1 + Sample12$hap2,
                      Sample13$hap1 + Sample13$hap2,
                      Sample14$hap1 + Sample14$hap2,
                      Sample15$hap1 + Sample15$hap2,
                      Sample18$hap1 + Sample18$hap2,
                      Sample19$hap1 + Sample19$hap2,
                      Sample20$hap1 + Sample20$hap2)
  
  GenoMatAfr = rbind(GetAncestrySpecCount(Sample9, 0), 
                    GetAncestrySpecCount(Sample10, 0), 
                    GetAncestrySpecCount(Sample11, 0), 
                    GetAncestrySpecCount(Sample12, 0), 
                    GetAncestrySpecCount(Sample13, 0), 
                    GetAncestrySpecCount(Sample14, 0), 
                    GetAncestrySpecCount(Sample15, 0), 
                    GetAncestrySpecCount(Sample18, 0), 
                    GetAncestrySpecCount(Sample19, 0), 
                    GetAncestrySpecCount(Sample20, 0)
  )
  
  
  GenoMatEur = rbind(GetAncestrySpecCount(Sample9, 1), 
                     GetAncestrySpecCount(Sample10, 1), 
                     GetAncestrySpecCount(Sample11, 1), 
                     GetAncestrySpecCount(Sample12, 1), 
                     GetAncestrySpecCount(Sample13, 1), 
                     GetAncestrySpecCount(Sample14, 1), 
                     GetAncestrySpecCount(Sample15, 1), 
                     GetAncestrySpecCount(Sample18, 1), 
                     GetAncestrySpecCount(Sample19, 1), 
                     GetAncestrySpecCount(Sample20, 1)
  )
  
  LAMat = rbind(Sample9$la1 + Sample9$la2,
                Sample10$la1 + Sample10$la2,
                Sample11$la1 + Sample11$la2,
                Sample12$la1 + Sample12$la2,
                Sample13$la1 + Sample13$la2,
                Sample14$la1 + Sample14$la2,
                Sample15$la1 + Sample15$la2,
                Sample18$la1 + Sample18$la2,
                Sample19$la1 + Sample19$la2,
                Sample20$la1 + Sample20$la2)
  
  Admprop = c( GetGlobAncestry(Sample9),
               GetGlobAncestry(Sample10),
               GetGlobAncestry(Sample11),
               GetGlobAncestry(Sample12),
               GetGlobAncestry(Sample13),
               GetGlobAncestry(Sample14),
               GetGlobAncestry(Sample15),
               GetGlobAncestry(Sample18),
               GetGlobAncestry(Sample19),
               GetGlobAncestry(Sample20))
  
  
  return(list(GenoMatTot = GenoMatTot, 
              GenoMatAfr = GenoMatAfr,
              GenoMatEur = GenoMatEur,
              LAMat = LAMat,
              Admprop = Admprop))
}





PlotPedigree <- function(){
  kindf <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),
                      dad =c(0, 0, 0, 1, 1, 0, 1, 0, 3,  3,  3,  5, 5,  7,  7,  7,  0,  16, 16, 16),
                      mom =c(0, 0, 0, 2, 2, 0, 2, 0, 4,  4,  4,  6, 6,  8,  8,  8,  0,  17, 17, 17),
                      sex =c(1, 0, 1, 0, 1, 0, 1, 0, 0,  1,  0,  1, 0,  1,  0,  1,  0,  0,  1,  0))
  return(plot(with(kindf, pedigree(id, mom, dad, sex))))
}



######################### Create a Genotype matrix ###############################

MakeGenotype <- function(nInd, nRel, MAFs, write = F, filename){
  
  famsize = 10 
  p = length(MAFs$maf1)
  

  if (nRel %% famsize != 0) stop('Family size is not a multiplication of n')
  
  # Generate Independent Individuals 
  GTot_Ind = matrix(data = NA, nrow = nInd, ncol = p)
  GAfr_Ind = matrix(data = NA, nrow = nInd, ncol = p)
  GEur_Ind = matrix(data = NA, nrow = nInd, ncol = p)
  LA_Ind   = matrix(data = NA, nrow = nInd, ncol = p)
  Admprop_Ind = c()
  
  for (i in 1:nInd){
    param1 = sample(c(8,12), size = 1)
    ind = GenerateAdm(MAFs$maf1, MAFs$maf2,  rbeta(1, param1, 20 - param1))
    GTot_Ind[i,] = ind$hap1 + ind$hap2
    GAfr_Ind[i,] = GetAncestrySpecCount(ind, 0)
    GEur_Ind[i,] = GetAncestrySpecCount(ind, 1)
    LA_Ind[i,]   = ind$la1 + ind$la2
    Admprop_Ind[i] = GetGlobAncestry(ind)
  }
  
  
  # Generate Family
  GTot_Rel = matrix(data = NA, nrow = nRel, ncol = p)
  GAfr_Rel = matrix(data = NA, nrow = nRel, ncol = p)
  GEur_Rel = matrix(data = NA, nrow = nRel, ncol = p)
  LA_Rel   = matrix(data = NA, nrow = nRel, ncol = p)
  Admprop_Rel = c()
  for (i in 1:(nRel/famsize)){
    param1 = sample(c(8,12), size = 1)
    fam = MakePedigree(MAFs, param1)
    GTot_Rel[(famsize * (i - 1) + 1):(famsize * i),] = fam$GenoMatTot
    GAfr_Rel[(famsize * (i - 1) + 1):(famsize * i),] = fam$GenoMatAfr
    GEur_Rel[(famsize * (i - 1) + 1):(famsize * i),] = fam$GenoMatEur
    LA_Rel[(famsize * (i - 1) + 1):(famsize * i),]   = fam$LAMat
    Admprop_Rel[(famsize * (i - 1) + 1):(famsize * i)] = fam$Admprop
  }

  
  GTot = rbind(GTot_Ind, GTot_Rel)
  GAfr = rbind(GAfr_Ind, GAfr_Rel)
  GEur = rbind(GEur_Ind, GEur_Rel)
  LA   = rbind(LA_Ind, LA_Rel)
  Admprop = c(Admprop_Ind, Admprop_Rel)
  
  # ensure the variant has variability and the allele frequency is between 0.1 - 0.9
  col2use = (apply(GTot, 2, var) != 0) & (apply(GAfr, 2, var) != 0) & (apply(GEur, 2, var) != 0) & (colMeans(GTot)/2 > 0.1) & (colMeans(GTot)/2 < 0.9)
  
  GTot = GTot[, col2use]
  GAfr = GAfr[, col2use]
  GEur = GEur[, col2use]
  LA   = LA[, col2use]
  
  row.names(GTot) = paste0("Sample", 1:nrow(GTot))
  row.names(GAfr) = paste0("Sample", 1:nrow(GAfr))
  row.names(GEur) = paste0("Sample", 1:nrow(GEur))
  row.names(LA) = paste0("Sample", 1:nrow(LA))
  
  
  if (write == T){
    WriteGenotype(GTot, paste0(filename, "_Tot.tsv")) 
    WriteGenotype(GAfr, paste0(filename, "_Afr.tsv"))
    WriteGenotype(GEur, paste0(filename, "_Eur.tsv"))
    WriteGenotype(LA, paste0(filename, "_LA.tsv"))
  }
  
  return(list(GTot = GTot,
              GAfr = GAfr, 
              GEur = GEur,
              LA   = LA,
              Admprop = Admprop))
  
}





MakeGRM <- function(nInd, nRel){
  famsize = 10
  if (nRel %% famsize != 0) stop('Family size is not a multiplication of n')
  GRM = matrix(data = 0, nrow = nInd + nRel, ncol = nInd + nRel)
  diag(GRM) = 1
  
  GRM_Rel = matrix(data = 0, nrow = nRel, ncol = nRel)
  kindf <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),
                      dad =c(0, 0, 0, 1, 1, 0, 1, 0, 3,  3,  3,  5, 5,  7,  7,  7,  0,  16, 16, 16),
                      mom =c(0, 0, 0, 2, 2, 0, 2, 0, 4,  4,  4,  6, 6,  8,  8,  8,  0,  17, 17, 17),
                      sex =c(1, 0, 1, 0, 1, 0, 1, 0, 0,  1,  0,  1, 0,  1,  0,  1,  0,  0,  1,  0))
  kinfam = kinship(with(kindf, pedigree(id, mom, dad, sex))) 
  kinfam = kinfam[c(9:15, 18:20), c(9:15, 18:20)] * 2
  for(i in 1:(nRel/famsize)){
    GRM_Rel[(famsize*(i-1)+1):(famsize*i), (famsize*(i-1)+1):(famsize*i)] = kinfam
  }
  
  GRM[(nInd + 1):(nInd + nRel), (nInd + 1):(nInd + nRel)] = GRM_Rel
  
  row.names(GRM) = paste0("Sample", 1:(nInd + nRel))
  colnames(GRM) = paste0("Sample", 1:(nInd + nRel))
  return(GRM)
}









