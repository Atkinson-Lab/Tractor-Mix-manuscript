#' Tractor score test
#'
#' @param obj The output form glmm.kin from GMMAT package
#' @param infiles a vector of input files in order (e.g. anc0.hapcount.txt, anc0.dosage.txt, anc1.dosage.txt); currently only support text file inputs.
#'
#' @return A list of summary statistics, based on users' input.
#' @export
TractorMix.score <- function(obj, infiles, outfiles) {
  # make sure infiles can take multiple input
  # the first line is header, will be used for sanity check
  # start for the 2nd line

  ################# step 1: read multiple files, and set up an iterator ##############
  # 1. use readLines to read all files
  # 2. save them in a list
  # 3. sanity check if the length is the same

  files = lapply(infiles,function(file){readLines(file)})
  filesName = lapply(infiles, function(file){
    fileSplit = unlist(strsplit(file, ".", fixed = T))
    return(paste(fileSplit[(length(fileSplit)-2):(length(fileSplit)-1)], collapse = "_"))
  })
  names(files) = filesName
  nSNP = length(files[[1]]) - 1


  ################# step 2: setup a result table ##############
  # 1. a list of empty dataframes, each one correspond to each ancestry in file name
  # 2. columns include: CHR, POS, ID, REF, ALT, P

  resDF = setNames(data.frame(matrix(data = NA, nrow = 0, ncol = (7 + 2 * length(infiles)))),
                                                       c("CHR", "POS", "ID", "REF", "ALT", "Chi2", "P", paste0("Eff_anc", 1:length(infiles)), paste0("Pval_anc", 1:length(infiles)) ))
  write.table(resDF, outfiles,  quote = F, row.names = F, sep = "\t")

  ################# step 3: calculate all necessary matricies ##############
  X = obj$X
  if(!is.null(obj$P)){
    ################################################################################################
    ################################### Full GRM ###################################################
    ################################################################################################


    W = diag(obj$fitted.values * (1 - obj$fitted.values))
    P = obj$P
    tildeT= diag(nrow(X)) - X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

    ################# step 4: iterate all files and calculate p values ##############
    for (i in 1:nSNP){
      # parsing a text file
      G = matrix(unlist(lapply(files, function(file){
        fileSplit = unlist(strsplit(file[i+1], "\t"))
        as.integer(fileSplit[6:length(fileSplit)])
      })), nrow = nrow(X), ncol = length(filesName))

      Gtilde = tildeT %*% G
      Score = t(G) %*% obj$scaled.residuals
      VarScore = t(Gtilde) %*% P %*% Gtilde
      chi2 = try(t(Score) %*%  solve(VarScore) %*% Score)

      if(class(chi2)[1] != "try-error") {
        p = pchisq(as.numeric(chi2), df = length(filesName), lower.tail = F)
        eff = t(Score) %*% solve(VarScore)
        se = sqrt(diag(solve(VarScore)))
        pval = sapply((eff/se)^2, function(teststats){pchisq(as.numeric(teststats), df = 1, lower.tail = F)})
        resDF[1, c("CHR", "POS", "ID", "REF", "ALT")] = unlist(strsplit(files[[1]][i+1], "\t"))[1:5]
        resDF[1, c("Chi2", "P")] = c(chi2, p)
        resDF[1, paste0("Eff_anc", 1:length(infiles))] = eff
        resDF[1, paste0("Pval_anc", 1:length(infiles))] = pval
        write.table(resDF,  outfiles, quote = F, row.names = F, col.names = F, append = T, sep = "\t")
      } else {
        resDF[1, c("CHR", "POS", "ID", "REF", "ALT")] = unlist(strsplit(files[[1]][i+1], "\t"))[1:5]
        resDF[1, c("Chi2", "P")] = c(NA, NA)
        resDF[1, paste0("Eff_anc", 1:length(infiles))] = rep(NA, length(infiles))
        resDF[1, paste0("Pval_anc", 1:length(infiles))] = rep(NA, length(infiles))
        write.table(resDF,  outfiles, quote = F, row.names = F, col.names = F, append = T, sep = "\t")
      }
    }
  } else {
    ################################################################################################
    ################################### Sparse GRM #################################################
    ################################################################################################


    obj$Sigma_iX = Matrix(obj$Sigma_iX, sparse = TRUE)  # this is solve(Sigma) %*% X
    obj$Sigma_i <- Matrix(obj$Sigma_i, sparse = TRUE)  # this is solve(Sigma)
    obj$cov = Matrix(obj$cov, sparse = TRUE)


    ################# step 4: iterate all files and calculate p values ##############
    for (i in 1:nSNP){
      # parsing a text file
      G = matrix(unlist(lapply(files, function(file){
        fileSplit = unlist(strsplit(file[i+1], "\t"))
        as.integer(fileSplit[6:length(fileSplit)])
      })), nrow = nrow(X), ncol = length(filesName))

      XSigma_iG = t(obj$Sigma_iX) %*% G
      VarScore = t(G) %*% obj$Sigma_i %*% G - t(XSigma_iG) %*% obj$cov %*% XSigma_iG
      Score = t(G) %*% obj$scaled.residuals
      chi2 = try(t(Score) %*%  solve(VarScore) %*% Score)

      if(class(chi2)[1] != "try-error") {
        p = pchisq(as.numeric(chi2), df = length(filesName), lower.tail = F)
        eff = t(Score) %*% solve(VarScore)
        se = sqrt(diag(solve(VarScore)))
        pval = sapply((eff/se)^2, function(teststats){pchisq(as.numeric(teststats), df = 1, lower.tail = F)})
        resDF[1, c("CHR", "POS", "ID", "REF", "ALT")] = unlist(strsplit(files[[1]][i+1], "\t"))[1:5]
        resDF[1, c("Chi2", "P")] = c(chi2, p)
        resDF[1, paste0("Eff_anc", 1:length(infiles))] = eff
        resDF[1, paste0("Pval_anc", 1:length(infiles))] = pval
        write.table(resDF,  outfiles, quote = F, row.names = F, col.names = F, append = T, sep = "\t")
      } else {
        resDF[1, c("CHR", "POS", "ID", "REF", "ALT")] = unlist(strsplit(files[[1]][i+1], "\t"))[1:5]
        resDF[1, c("Chi2", "P")] = c(NA, NA)
        resDF[1, paste0("Eff_anc", 1:length(infiles))] = rep(NA, length(infiles))
        resDF[1, paste0("Pval_anc", 1:length(infiles))] = rep(NA, length(infiles))
        write.table(resDF,  outfiles, quote = F, row.names = F, col.names = F, append = T, sep = "\t")
      }
    }

  }


}

