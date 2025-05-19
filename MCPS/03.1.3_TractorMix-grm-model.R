# Define the custom paste function
"%&%" <- function(a, b) paste0(a, b)
library(GENESIS)
library(GMMAT)
library(GWASTools)
library(SNPRelate)
library(dplyr)
library(data.table)
library(zoo)
library(optparse)
library(cli)

proj_dir <- "/well/emberson/users/pom143/projects/LAI_TractorMix/"
out_dir <- proj_dir %&% "output-files/"
bmi_df <- fread(proj_dir %&% "input-files/extracted-phenotypes-all-2.5perc-bmi_renamed.txt", header = T)

spinny <- make_spinner(
    which = "line",
    template = "{spin} Running process."
)

make_grms <- function(analysis_name){
  ## analysis_name is either "analysis-set-1_bmi" or "analysis-set-2_hypertension"
  cli::cli_alert("Running steps for analysis: " %&% analysis_name)
  if (analysis_name == "analysis-set-1_bmi"){
    pheno_df <- bmi_df
  } else{
    stop("Analysis name must be 'analysis-set-1_bmi'")
  }

  load(out_dir %&% analysis_name %&% "_pcrelate_results_33500.RData")
  cli::cli_alert("Creating GRM matrices...")
  cli::cli_alert("Full matrix...")
  iid <- pheno_df$IID
  spinny$spin()
  GRM <- pcrelateToMatrix(pcrelate_r, sample.include = iid, scaleKin = 2)[iid, iid]
  spinny$finish()
  cli::cli_alert("Saving full matrix")
  save(GRM, file = out_dir %&% analysis_name %&% "_GRM_c.RData")

  cli::cli_alert("Sparse matrix...")
  spinny$spin()
  GRM_sparse <- GRM
  GRM_sparse[GRM_sparse < 0.05] <- 0
  GRM_sparse <- as(GRM_sparse, "sparseMatrix")
  spinny$finish()
  cli::cli_alert("Saving sparse matrix")
  save(GRM_sparse, file = out_dir %&% analysis_name %&% "_GRM_sparse_c.RData")
  cli::cli_alert_success("GRM process complete.")
  return(list(GRM, GRM_sparse))
}

grm_list_1 <- make_grms("analysis-set-1_bmi")


# Fit the null model with GMMAT
# 'phenoDf' dataframe with ID, Pheno, age, sex, and PCs from PC-AiR

# Continuous phenotype
Model_Null_continuous_gm <- glmmkin(fixed = BMI ~ SEX + AGE + Score_AMR, data = bmi_df, 
    id = "IID", kins = grm_list_1[[1]], family = gaussian(link = "identity"))
Model_Null_continuous_gm_sp <- glmmkin(fixed = BMI ~ SEX + AGE + Score_AMR, data = bmi_df, 
    id = "IID", kins = grm_list_1[[2]], family = gaussian(link = "identity"))
save(Model_Null_continuous_gm, file = out_dir %&% "model_null_bmi_full-matrix-33500.RData")
save(Model_Null_continuous_gm_sp, file = out_dir %&% "model_null_bmi_sparse-matrix-33500.RData")



