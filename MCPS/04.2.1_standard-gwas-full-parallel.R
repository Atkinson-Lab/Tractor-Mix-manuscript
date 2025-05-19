# Define the custom paste function
"%&%" <- function(a, b) paste0(a, b)
library(GENESIS)
library(GMMAT)
library(GWASTools)
library(dplyr)
library(data.table)
library(zoo)
library(cli)
library(optparse)
library(BiocParallel)

# Define parallel parameters
bpparam <- BiocParallel::SnowParam(workers = 4, type = "SOCK")

# Define command-line options
option_list = list(
  make_option(c("--chr"), type="integer", default=NULL, help="Chromosome number")
)

# Parse command-line options
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check if --chr option is provided
if (is.null(opt$chr)) {
  stop("Please provide the chromosome number using the --chr option.")
}

chr <- opt$chr

# Function to run GWAS for BMI or Hypertension
run_gwas <- function(analysis_name, chr) {
    "%&%" <- function(a, b) paste0(a, b)
    cli::cli_alert("Running steps for analysis: " %&% analysis_name)

    proj_dir <- "/well/emberson/users/pom143/projects/LAI_TractorMix/"
    inp_dir <- proj_dir %&% "output-files/"
    bed_dir <- "/well/emberson/users/pom143/projects/LAI_TractorMix/output-files/chr_autosomes_sample/"
    out_dir <- "/well/emberson/users/pom143/projects/LAI_TractorMix/output-files/standard_gwas/"

    # Load the appropriate models based on analysis_name
    if (analysis_name == "analysis-set-1_bmi") {
        load(paste0(inp_dir, "model_null_bmi_full-matrix-33500.RData"))
        model_null <- Model_Null_continuous_gm  # Ensure this is the correct object name
        result_file <- out_dir %&% "bmi_chr" %&% chr %&% "_gwas_results_full_33500.txt"
    } else {
        stop("Invalid analysis name")
    }

    # Generate the file prefix for the chromosome-specific bed file
    geno_prefix <- bed_dir %&% "merged_hgdp-1kg-mcps_autosomes_33500_chr" %&% chr

    # Check if the required files exist before proceeding
    if (!file.exists(geno_prefix %&% ".bed")) {
        stop("Error: The input BED file does not exist at: " %&% geno_prefix %&% ".bed")
    }
    if (!file.exists(geno_prefix %&% ".bim")) {
        stop("Error: The input BIM file does not exist at: " %&% geno_prefix %&% ".bim")
    }
    if (!file.exists(geno_prefix %&% ".fam")) {
        stop("Error: The input FAM file does not exist at: " %&% geno_prefix %&% ".fam")
    }

    cli::cli_alert("Running GWAS Score Test for Chromosome " %&% chr %&% " on " %&% analysis_name %&% " analysis...")

    # Run the GWAS using the loaded model and chromosome-specific bed file
    GMMAT::glmm.score(model_null, infile = geno_prefix, outfile = result_file, missing.method = "omit", center = F)
    
    cli::cli_alert_success("Test complete for Chromosome " %&% chr %&% ". Results saved to: " %&% result_file)
}

# List of analyses to loop through
analyses <- c("analysis-set-1_bmi")

# Use bplapply to parallelize the loop across analyses and chromosomes
results <- BiocParallel::bplapply(analyses, function(analysis_name) run_gwas(analysis_name, chr), BPPARAM = bpparam)

