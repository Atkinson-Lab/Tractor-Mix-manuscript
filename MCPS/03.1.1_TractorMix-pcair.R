# Define the custom paste function
"%&%" <- function(a, b) paste0(a, b)
library(GENESIS)
## library(GMMAT)
library(GWASTools)
library(SNPRelate)
library(dplyr)
library(data.table)
library(zoo)
library(cli)

proj_dir <- "/well/emberson/users/pom143/projects/LAI_TractorMix/"
out_dir <- proj_dir %&% "output-files/"

spinny <- make_spinner(
    which = "line",
    template = "{spin} Running process."
)

run_pcair <- function(analysis_name){
  ## analysis_name is either "analysis-set-1_bmi" or "analysis-set-2_hypertension"
  cli::cli_alert("Running steps for analysis: " %&% analysis_name)
  bed_file <- out_dir %&% analysis_name %&% ".bed"
  bim_file <- out_dir %&% analysis_name %&% ".bim"
  fam_file <- out_dir %&% analysis_name %&% ".fam"
  gds_file <- out_dir %&% analysis_name %&% ".gds"
  kin0_file <- out_dir %&% "king/" %&% analysis_name %&% ".kin0"
  kin_file <- out_dir %&% "king/" %&% analysis_name %&% ".kin"

  cli::cli_alert("Convert BED files to GDS format...")
  spinny$spin()
  snpgdsBED2GDS(bed.fn = bed_file, bim.fn = bim_file, fam.fn = fam_file, out.gdsfn = gds_file)
  gds <- snpgdsOpen(gds_file)
  spinny$finish()
  cli::cli_alert_success("Conversion complete.")

  cli::cli_alert("LD pruning GDS file...")
  spinny$spin()
  snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
    ld.threshold=sqrt(0.1), verbose=FALSE)
  spinny$finish()
  cli::cli_alert_success("Pruning complete.")

  cli::cli_alert("Saving data...")
  spinny$spin()
  save(snpset, file = out_dir %&% analysis_name %&% "_snpset.RData")
  spinny$finish()
  cli::cli_alert_success("Save complete.")
  pruned <- unlist(snpset, use.names=FALSE)

  cli::cli_alert("Creating KING matrix...")

  sample_names <- read.gdsn(index.gdsn(gds, "sample.id"))
  
  # Load kin0 file using fread for memory efficiency
  kinship_data <- fread(kin0_file)
  
  # Combine ID1 and ID2 to get all unique individuals
  all_ids <- unique(c(kinship_data$ID1, kinship_data$ID2))

  # Ensure that the IDs in the kinship data match the sample names from the GDS
  if (!all(all_ids %in% sample_names)) {
    stop("Some IDs in the kinship file do not match the sample names in the GDS file.")
  }

  id_map <- setNames(1:length(sample_names), sample_names)  # Create a map for ID to index

  # Convert IDs to numeric indices for sparse matrix representation
  id1 <- id_map[as.character(kinship_data$ID1)]
  id2 <- id_map[as.character(kinship_data$ID2)]
  
  # Create sparse matrix representation of the kinship data
  sparse_kinship <- sparseMatrix(
    i = id1,  # Row indices
    j = id2,  # Column indices
    x = kinship_data$Kinship,  # Kinship values
    dims = c(length(sample_names), length(sample_names)),  # Make it square
    symmetric = TRUE  # Kinship matrix is symmetric
  )

  # Assign row and column names to match sample names from the GDS
  rownames(sparse_kinship) <- sample_names
  colnames(sparse_kinship) <- sample_names

  KINGmat2 <- as(sparse_kinship, "sparseMatrix")
  cli::cli_alert_success("Kinship matrix loaded as sparse matrix.")

  cli::cli_alert("Memory usage before PC-Air:")
  print(gc())
  
  cli::cli_alert("Running PC-Air...")
  spinny$spin()
  tryCatch({
    pcair_res <- pcair(gds, kinobj = KINGmat2, divobj = KINGmat2, snp.include = pruned, num.cores = 16)
    save(pcair_res, file = out_dir %&% analysis_name %&% "_pcair_results.RData")
    spinny$finish()
    cli::cli_alert_success("PC-AIR complete.")
  }, error = function(e) {
    spinny$finish()
    cli::cli_alert_danger("Error during PC-AIR: " %&% e$message)
    snpgdsClose(gds)
    stop("PC-AIR failed.")
  })

  snpgdsClose(gds)
  rm(gds, KINGmat2, pcair_res, pruned)
  gc()
}

run_pcair("analysis-set-1_bmi")


