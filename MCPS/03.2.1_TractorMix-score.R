# Define the custom paste function
"%&%" <- function(a, b) paste0(a, b)

# Load required libraries
library(GENESIS)
library(GMMAT)
library(GWASTools)
library(SNPRelate)
library(dplyr)
library(data.table)
library(zoo)
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("--chr"), type = "integer", default = NULL, help = "Chromosome number")
)

# Parse command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if --chr option is provided
if (is.null(opt$chr)) {
  stop("Please provide the chromosome number using the --chr option.")
}

chr <- opt$chr

# Define file paths
base_dir <- "/well/emberson/users/pom143/projects/LAI_TractorMix/"
file_rfmix <- "/well/emberson/users/pom143/projects/LAI_TractorMix/output-files/dosages/"
output_file <- "/well/emberson/users/pom143/projects/LAI_TractorMix/output-files/tractor-mix/"
tracmix_script <- "/well/emberson/users/pom143/projects/LAI_TractorMix/src/development/"

# Source external script
source(paste0(tracmix_script, "TractorMix.score.update2.R"))

# Function to load RData and return specific objects
load_model_data <- function(full_matrix_file, sparse_matrix_file) {
  # Load the full matrix
  load(full_matrix_file)
  full_model_name <- ls()[grep("Model_Null_", ls())]  # Find the correct full model object
  
  # Load the sparse matrix
  load(sparse_matrix_file)
  sparse_model_name <- ls()[grep("Model_Null_", ls())]  # Find the correct sparse model object
  
  # Print loaded objects for debugging
  print(paste("Full model loaded:", full_model_name))
  print(paste("Sparse model loaded:", sparse_model_name))
  
  # Return both models
  return(list(
    full_model = get(full_model_name[1]),
    sparse_model = get(sparse_model_name[1])
  ))
}

# Function to run TractorMix.score
run_tractormix <- function(model, phenotype, grm_type) {
  infiles <- c(
    paste0(file_rfmix, "mcps_sample_amr_eur_phased_", phenotype, "_chr", chr, ".anc0.dosage.txt"),
    paste0(file_rfmix, "mcps_sample_amr_eur_phased_", phenotype, "_chr", chr, ".anc1.dosage.txt")
  )
  
  outfiles <- paste0(output_file, "tractorMix_", phenotype, "_", grm_type, "-grm_chr", chr, ".txt")
  
  TractorMix.score(
    obj = model,  
    infiles = infiles,
    outfiles = outfiles
  )
}

# Load BMI models (ensure proper file paths)
bmi_models <- load_model_data(
  full_matrix_file = paste0(base_dir, "output-files/model_null_bmi_full-matrix-33500.RData"),
  sparse_matrix_file = paste0(base_dir, "output-files/model_null_bmi_sparse-matrix-33500.RData")
)

# Check if the models are loaded correctly
if (is.null(bmi_models$full_model)) {
  stop("BMI full model not loaded correctly.")
}
if (is.null(bmi_models$sparse_model)) {
  stop("BMI sparse model not loaded correctly.")
}

# Run TractorMix for BMI (full and sparse)
run_tractormix(bmi_models$full_model, "bmi", "full")
run_tractormix(bmi_models$sparse_model, "bmi", "sparse")





