#' Tractor score test
#'
#' @param obj The output form glmm.kin from GMMAT package
#' @param infiles a vector of input file names in order (e.g. anc0.dosage.txt, anc1.dosage.txt)
#' @return A list of summary statistics, based on users' input.
#' @export
#' 

# version: 0.0.1
script_version <- "0.0.1"
cat(paste("TractorMix.score Script Version :", script_version), "\n")



TractorMix.score <- function(obj, infiles, outfiles, n_core = 4, chunk_size = 2048) {
  
  ## install required packages
  required_packages <- c("Matrix", "data.table", "doParallel", "foreach", "dplyr", "abind")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      suppressMessages(library(pkg, character.only = TRUE))
    }
  }
  
  # determine if n_core is correct
  if(n_core>detectCores()) {
    stop(paste0("Error: n_core must be less than ", detectCores(),"!"))
  }
  if (n_core > 1) {
    n_core = n_core - 1
  } else if (n_core <= 0) {
    stop(paste0("Error: --nthreads must be 1 or more!"))
  }
  
  # number of ancestries
  n_anc = length(infiles)
  

  # use wc -l to find the number of SNPs in each file
  # ASSUME: dosage files are equal length
  wc_op <- paste("wc -l", infiles[1]) %>%
    system(intern = T) %>%
    trimws("left") %>% 
    strsplit(" ")
  
  n_SNPs = wc_op[[1]][1] %>% as.numeric() -1
  
  
  # setup a result table 
  # 2-way CHR, POS, ID, REF, ALT, Chi2, P, Eff_anc0, Eff_anc1, Pval_anc0, Pval_anc1, AC_anc0, AC_anc1
  # AC might explain FP inflation for binary outcome
  resDF = setNames(data.frame(matrix(data = NA, nrow = 0, ncol = (7 + 4 * n_anc))),
                   c("CHR", "POS", "ID", "REF", "ALT", "Chi2", "P", 
                     paste0("Eff_anc", 0:(n_anc-1)), 
                     paste0("SE_anc", 0:(n_anc-1)), 
                     paste0("Pval_anc", 0:(n_anc-1)),
                     paste0("AC_anc", 0:(n_anc-1))))
  write.table(resDF, outfiles,  quote = F, row.names = F, sep = "\t")
  
  
  # retrieve file names, and check if consistent
  files_colnames = lapply(infiles, function(infile){colnames(fread(input = infile, sep = "\t", skip = 0, nrows = 1))})
  if (Reduce(function(x, y) identical(x, y), files_colnames)){
    sample_id = files_colnames[[1]][6:length(files_colnames[[1]])]
  }else{
    stop("Sample order in dosage files do not match!")
  }
  # ASSUME: obj$id_include is a subset of samples in dosage files
  sample2Use = match(obj$id_include, sample_id)
  if (any(is.na(sample2Use))){
    stop("The null model contains samples that are not in the dosage files; 
         remove unmatched samples before fitting the null model")
  }
  
  # define column class
  colclasses = c(rep("character", 5), rep("integer", length(sample_id)))
  
  
  if(!is.null(obj$P)){
    ################################################################################################
    ################################### Full GRM ###################################################
    ################################################################################################
    # Score test; shared across all SNPs
    X = obj$X
    W = diag(obj$fitted.values * (1 - obj$fitted.values))
    P = obj$P
    tildeT= diag(nrow(X)) - X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    
    
    # Parallel
    cl <- makeCluster(n_core)
    registerDoParallel(cl)
    pb <- txtProgressBar(min = 0, max = ceiling(n_SNPs / chunk_size), style = 3)
    
    chunk_idx = 0
    while (chunk_idx < ceiling(n_SNPs / chunk_size)){
      files = lapply(infiles, 
                     function(infile){fread(input = infile, 
                                            sep = "\t", 
                                            skip = 1 + chunk_idx * chunk_size, 
                                            nrows = chunk_size, 
                                            colClasses = colclasses)})
      
      meta_info = files[[1]][,1:5]
      geno_info = abind(lapply(files, function(file){ unname(as.matrix(file[, 6:(length(sample_id) + 5)])[,sample2Use]) }), along = 3)
      
      
      # parallelization within each chunk
      n_SNPs_chunk = nrow(files[[1]])
      sumstats_chunk = foreach(i = 1:n_SNPs_chunk, .combine = rbind, .inorder = TRUE) %dopar% {
        # parse genotype
        G = geno_info[i,,]
        Score = t(G) %*% obj$scaled.residuals
        VarScore = t(G) %*% P %*% G
        joint_chi2 = try(t(Score) %*%  solve(VarScore) %*% Score, silent=TRUE)
        if(class(joint_chi2)[1] != "try-error") { # VarScore is invertable
          joint_pval = pchisq(as.numeric(joint_chi2), df = n_anc, lower.tail = F)
          anc_eff = t(Score) %*% solve(VarScore)
          anc_se = sqrt(diag(solve(VarScore)))
          anc_pval = sapply((anc_eff/anc_se)^2, function(teststats){pchisq(as.numeric(teststats), df = 1, lower.tail = F)})
          res = c(round(joint_chi2, 5), 
                  signif(joint_pval, 5), 
                  round(anc_eff, 5), 
                  round(anc_se, 5), 
                  signif(anc_pval, 5), 
                  colSums(G))
        } else { # VarScore is not invertable
          res = c(rep(NA, 2 + 3 * n_anc), 
                  colSums(G))
        }
        names(res) = NULL
        return(res)
      }
      
      write.table(cbind(meta_info, sumstats_chunk), outfiles,  quote = F, row.names = F, sep = "\t", append = T, col.names = F)
      chunk_idx = chunk_idx + 1
      setTxtProgressBar(pb, chunk_idx) 
    }
    close(pb)
    stopCluster(cl)
    
  } else {
    ################################################################################################
    ################################### Sparse GRM #################################################
    ################################################################################################
    obj$Sigma_iX = Matrix(obj$Sigma_iX, sparse = TRUE)  # this is solve(Sigma) %*% X
    obj$Sigma_i <- Matrix(obj$Sigma_i, sparse = TRUE)  # this is solve(Sigma)
    obj$cov = Matrix(obj$cov, sparse = TRUE)
    
    # Parallel
    cl <- makeCluster(n_core)
    registerDoParallel(cl)
    pb <- txtProgressBar(min = 0, max = ceiling(n_SNPs / chunk_size), style = 3)
    
    chunk_idx = 0
    while (chunk_idx < ceiling(n_SNPs / chunk_size)){
      files = lapply(infiles, 
                     function(infile){fread(input = infile, sep = "\t", 
                                            skip = 1 + chunk_idx * chunk_size, 
                                            nrows = chunk_size, 
                                            colClasses = colclasses)})
      
      meta_info = files[[1]][,1:5]
      geno_info = abind(lapply(files, function(file){ unname(as.matrix(file[, 6:(length(sample_id) + 5)])[,sample2Use]) }), along = 3)
      
      # parallelization within each chunk
      n_SNPs_chunk = nrow(files[[1]])
      sumstats_chunk = foreach(i = 1:n_SNPs_chunk, .combine = rbind, .inorder = TRUE, .packages=c("Matrix") ) %dopar% {
        # parse genotype
        G = geno_info[i,,]
        XSigma_iG = t(obj$Sigma_iX) %*% G
        VarScore = t(G) %*% obj$Sigma_i %*% G - t(XSigma_iG) %*% obj$cov %*% XSigma_iG
        Score = t(G) %*% obj$scaled.residuals
        joint_chi2 = try(t(Score) %*%  solve(VarScore) %*% Score, silent = T)

        if(class(joint_chi2)[1] != "try-error") { # VarScore is invertable
          joint_pval = pchisq(as.numeric(joint_chi2), df = n_anc, lower.tail = F)
          anc_eff = t(Score) %*% solve(VarScore)
          anc_se = sqrt(diag(solve(VarScore)))
          anc_pval = sapply((anc_eff/anc_se)^2, function(teststats){pchisq(as.numeric(teststats), df = 1, lower.tail = F)})
          res = c(round(as.numeric(joint_chi2), 5),
                  signif(as.numeric(joint_pval), 5),
                  round(as.numeric(anc_eff), 5),
                  round(as.numeric(anc_se), 5),
                  signif(as.numeric(anc_pval), 5),
                  colSums(G))
        } else { # VarScore is not invertable
          res = c(rep(NA, 2 + 3 * n_anc),
                  colSums(G))
        }
        names(res) = NULL
        return(res)
      }
      
      write.table(cbind(meta_info, sumstats_chunk), outfiles,  quote = F, row.names = F, sep = "\t", append = T, col.names = F)
      chunk_idx = chunk_idx + 1
      setTxtProgressBar(pb, chunk_idx) 
      }
    close(pb)
    stopCluster(cl)
    }
    
  }
  
  
  
  
