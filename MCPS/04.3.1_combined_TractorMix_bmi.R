"%&%" <- function(a,b) paste0(a,b)
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("qqman"))

# Set the output directory
inp.dir <- "/well/emberson/users/pom143/projects/LAI_TractorMix/output-files/tractor-mix/"
plot.dir <- "/well/emberson/users/pom143/projects/LAI_TractorMix/output-files/tractor-mix/plots/"
out.dir <- "/well/emberson/users/pom143/projects/LAI_TractorMix/output-files/tractor-mix/"

# Check if directories exist, create plot directory if not
if (!dir.exists(inp.dir)) {
  stop("The specified input directory does not exist: " %&% inp.dir)
}

if (!dir.exists(plot.dir)) {
  dir.create(plot.dir)
}

if (!dir.exists(out.dir)) {
  stop("The specified output directory does not exist: " %&% out.dir)
}

# List of base file names to process
file_bases <- c("tractorMix_bmi_full-grm_chr")

# List of columns to loop through
columns_to_process <- c("P", "Pval_anc0", "Pval_anc1")

# Function to merge all chromosomes
merge_chromosomes <- function(file_base) {
  merged_data <- NULL
  for (chr in 1:22) {
    file.name <- file_base %&% chr %&% ".txt"
    file.path <- inp.dir %&% file.name
    
    if (!file.exists(file.path)) {
      stop("The specified file does not exist: " %&% file.path)
    }
    
    chr_data <- fread(file.path)
    merged_data <- rbind(merged_data, chr_data)
  }
  return(merged_data)
}

# Loop over each file base
for (file_base in file_bases) {
  print("Merging chromosomes for: " %&% file_base)
  merged_data <- merge_chromosomes(file_base)
  
  for (col in columns_to_process) {
    print("Processing column: " %&% col)
    
    if (col == "P") {
      merged_data$Z <- sqrt(merged_data$Chi2)
      inflation.factor <- median((na.omit(merged_data$Z)^2) / qchisq(1/2, df=2))
      output.file <- out.dir %&% file_base %&% "_combined_P.txt"
      
      # Select relevant columns for output
      gwas_sub.df <- merged_data %>%
        select(CHR, POS, ID, REF, ALT, Chi2, Z, P)
      
      # Rename columns
      colnames(gwas_sub.df) <- c("CHR", "POS", "ID", "REF", "ALT", "Chi2", "Z", "P")
      
    } else {
      anc <- gsub("Pval_", "", col)  # Extract anc0 or anc1
      eff_col <- "Eff_" %&% anc
      se_col <- "SE_" %&% anc
      
      merged_data$Z <- merged_data[[eff_col]] / merged_data[[se_col]]
      inflation.factor <- median((na.omit(merged_data$Z)^2) / qchisq(1/2, df=1))
      output.file <- out.dir %&% file_base %&% "_combined_" %&% anc %&% ".txt"
      
      # Select relevant columns for output
      gwas_sub.df <- merged_data %>%
        select(CHR, POS, ID, REF, ALT, Z, !!sym(col), !!sym(eff_col), !!sym(se_col))
      
      # Rename columns to reflect ancestry-specific effects
      colnames(gwas_sub.df) <- c("CHR", "POS", "ID", "REF", "ALT", "Z", col, eff_col, se_col)
    }
    
    print("Writing output file: " %&% output.file)
    write.table(gwas_sub.df, file=output.file, sep="\t", col.names=T, row.names=F, quote=F)
    
    ## Create subsetted data frame to use for generating Manhattan plot
    print("Creating subsetted data frame to be used for Manhattan plot...")
    sig <- 5e-8  # Genome-wide significance threshold
    sig.df <- filter(gwas_sub.df, !!sym(col) < sig)
    
    # If no significant variants, continue to next column
    if (nrow(sig.df) == 0) {
      print("No significant variants for column: " %&% col)
      next
    }
    
    # Sample non-significant variants for plotting
    to_sample <- 1e6 - nrow(sig.df)
    null_count <- (to_sample / 2) %>% round(.)
    
    notsig.df <- filter(gwas_sub.df, !(ID %in% sig.df$ID))
    notsig.dfa <- filter(notsig.df, !!sym(col) >= 0.001, !!sym(col) < 0.05)
    notsig.dfb <- filter(notsig.df, !!sym(col) > 0.05)
    
    # Ensure sample size does not exceed the population size
    keep_index1 <- sample(1:nrow(notsig.dfa), size=min(null_count, nrow(notsig.dfa)), replace = F)
    keep_index2 <- sample(1:nrow(notsig.dfb), size=min(null_count, nrow(notsig.dfb)), replace = F)
    
    null_df1 <- notsig.dfa[keep_index1, ]
    null_df2 <- notsig.dfb[keep_index2, ]
    
    plot.df <- rbind(null_df1, null_df2, sig.df)
    manhattan.input.file <- out.dir %&% "manhattan-input_" %&% file_base %&% "_" %&% col %&% ".txt"
    print("Writing output file: " %&% manhattan.input.file)
    write.table(plot.df, file=manhattan.input.file, sep="\t", col.names=T, row.names=F, quote=F)
    
    ## Report significant variant counts and inflation factor
    print("There are " %&% dim(sig.df)[1] %&% " significant variants.\n")
    print("Observed lambda inflation factor: " %&% inflation.factor %&% "\n")
    
    # Generate QQ-plot
    qq.plot.file <- plot.dir %&% "qq_" %&% file_base %&% "_" %&% col %&% ".jpeg"
    print("Saving QQ-plot image file: " %&% qq.plot.file)
    png(qq.plot.file, type = "cairo")
    qq(merged_data[[col]], main="Inflation Factor: " %&% round(inflation.factor, digits=3))
    dev.off()
    
    if (dim(sig.df)[1] > 0) {
      genome.wide.sig.file <- out.dir %&% "genome-wide-significant_" %&% file_base %&% "_" %&% col %&% ".txt"
      print("Printing Genome-wide significant associations")
      write.table(filter(sig.df, !!sym(col) < sig), file=genome.wide.sig.file, sep="\t", col.names=T, row.names=F, quote=F)
    }
  }
}

