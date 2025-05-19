## Developer: Alejandra Vergara-Lope & Jason Matthew Torres
## module load R/4.3.2-gfbf-2023a
## module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

set.seed(1)

# Define a string concatenation function
"%&%" <- function(a, b) paste0(a, b)

# Load necessary libraries
library(data.table)
library(dplyr)

# Define directories
base_dir <- "/well/emberson/users/bjk420/projects/"
work_dir <- base_dir %&% "LAI_TractorMix/"
input_dir <- work_dir %&% "input-files/"
ancestry_dir <- base_dir %&% "ancestry-proportions/output-files/"
output_dir <- "/well/emberson/users/pom143/projects/LAI_TractorMix/input-files/"

# Load ancestry and phenotype data
anc_df <- fread(ancestry_dir %&% "ancestry-proportions-k4.txt", header = TRUE)
pheno_df <- fread(input_dir %&% "extracted-phenotypes.txt", header = TRUE)
# Function to filter and write output using write.table
filter_and_save <- function(pheno_data, suffix) {
  bmi_notmiss_df <- filter(pheno_data, !is.na(BMI))  # Filter out individuals with missing BMI
  print(dim(bmi_notmiss_df))
  write.table(bmi_notmiss_df, file = output_dir %&% suffix %&% ".txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Combine phenotypes and ancestry data for BMI analysis
suffix1 <- "extracted-phenotypes-all-bmi"
pheno_filt <- filter(pheno_df, BASE_DIABETES == 0, BASE_HBA1C <= 6)
pheno_sub_df <- left_join(pheno_filt, anc_df, by = "IID")

# Filter for non-missing BMI and save
filter_and_save(pheno_sub_df, suffix1)
print(dim(pheno_sub_df))

# Apply ancestry thresholds for EAS and AFR ancestry
process_ancestry_thresholds <- function(anc_df, pheno_df, thresholds, suffixes) {
  pheno_filt_df <- filter(pheno_df, BASE_DIABETES == 0, BASE_HBA1C <= 6)
  
  # Iterate over the thresholds and suffixes to process each case
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    suffix <- suffixes[i]
    
    # Apply ancestry thresholds and merge
    anc_df$Score_eas_afr <- anc_df$Score_EAS + anc_df$Score_AFR
    anc_filt_df <- filter(anc_df, Score_eas_afr <= threshold, population == "MCPS")
    
    pheno_sub_df <- left_join(filter(pheno_filt_df, IID %in% anc_filt_df$IID), anc_filt_df, by = "IID")
    
    # Filter for non-missing BMI and save
    filter_and_save(pheno_sub_df, suffix)
    print(dim(pheno_sub_df))
  }
}

# Define thresholds and suffixes for ancestry filtering
thresholds <- c(0.02, 0.025, 0.05)
suffixes <- c("extracted-phenotypes-all-2perc-bmi", "extracted-phenotypes-all-2.5perc-bmi", "extracted-phenotypes-all-5perc-bmi")

# Process ancestry thresholds
process_ancestry_thresholds(anc_df, pheno_df, thresholds, suffixes)
