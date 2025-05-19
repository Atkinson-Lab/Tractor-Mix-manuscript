"%&%" <- function(a,b) paste0(a,b)
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tools"))

# Functions
process_gwas_df <- function(gwas.dat){
  gwas.dat <- as.data.frame(gwas.dat)
  chrom.levels <- c(1:22)
  gwas.dat$CHR <- factor(gwas.dat$CHR, levels = chrom.levels)
  gwas.dat <- arrange(gwas.dat, CHR)
  nCHR <- length(unique(gwas.dat$CHR))
  gwas.dat$POScum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(gwas.dat$CHR)){
    nbp[i] <- max(gwas.dat[gwas.dat$CHR == i,]$POS)
    gwas.dat[gwas.dat$CHR == i, "POScum"] <- gwas.dat[gwas.dat$CHR == i, "POS"] + s
    s <- s + nbp[i]
  }
  return(gwas.dat)
}

annotate_genes_to_sig_snps <- function(gwas.dat, sig = 5e-8, gene.gr){
  sig.df <- filter(gwas.dat, P <= sig)
  nearest.gene <- c()
  for (i in 1:nrow(sig.df)){
    sig.gr <- GRanges(seqnames = sig.df$CHR[i],
                      ranges = IRanges(start = sig.df$POS[i], end = sig.df$POS[i]))
    ng <- names(gene.gr)[nearest(sig.gr, gene_gr)]
    nearest.gene <- append(nearest.gene, ng)
  }
  sig.df$nearest.gene <- nearest.gene
  out.df <- data.frame()
  for (gene in unique(sig.df$nearest.gene)){
    sub.df <- filter(sig.df, nearest.gene == gene) %>% arrange(P)
    out.df <- rbind(out.df, sub.df[1,])
  }
  return(out.df)
}

manhat <- function(gwas.dat, p_col, sig = 5e-8, gene.gr){
  nCHR <- length(unique(gwas.dat$CHR))
  axis.set <- suppressMessages(gwas.dat %>% group_by(CHR) %>%
                                 summarize(center = (max(POScum) + min(POScum)) / 2))
  ylim <- abs(floor(log10(min(na.omit(gwas.dat[[p_col]]))))) + 2
  gwas.dat$P <- gwas.dat[[p_col]]
  sig.df <- annotate_genes_to_sig_snps(gwas.dat, sig, gene.gr)
  sig.df <- arrange(sig.df, P)
  plt <- ggplot(gwas.dat, aes(x = POScum, y = -log10(P))) +
    geom_point(aes(color = as.factor(CHR)), alpha = 0.75, size = 1.25) +
    geom_hline(yintercept = -log10(5e-8), color = viridis(10)[1], linetype = "dashed") +
    scale_color_manual(values = rep(c(color_1, color_2), nCHR)) +
    scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
    labs(x = NULL, y = "-log10(p)") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0, size = 8, vjust = 0.5))
  
  if (nrow(sig.df) > 0){
    if (nrow(sig.df) <= 100){
      full.plt <- plt +
        geom_text_repel(data = sig.df,
                        aes(x = POScum, y = -log10(P), label = nearest.gene, size = -log10(P))) +
        scale_size(range = c(1.5, 4))
    } else {
      full.plt <- plt +
        geom_text_repel(data = sig.df[1:100,],
                        aes(x = POScum, y = -log10(P), label = nearest.gene, size = -log10(P))) +
        scale_size(range = c(1.5, 4))
    }
  } else {
    full.plt <- plt
  }
  return(list(full.plt, plt, sig.df, gwas.dat))
}

# Set up directories
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Please provide the input and output directories as arguments.")
}
input_dir <- args[1]
output_dir <- args[2]
output_dir <- ifelse(str_sub(output_dir,-1)=="/",output_dir,output_dir%&%"/")
plot_dir <- output_dir %&% "plots/"
if (!dir.exists(plot_dir)){ dir.create(plot_dir) }

if (length(args)==4){
  color_1 <- args[3]
  color_2 <- args[4]
} else{
  color_1 <- "#31688EFF"
  color_2 <- "#1F9E89FF"
}

# Load the gene annotations
message("Reading in ENSEMBL GRCh38 gene annotations...")
ens_dir <- "/well/emberson/shared/reference_datasets/ensembl-genes_build-38/"
ens_df <- import(ens_dir %&% "Homo_sapiens.GRCh38.105.gtf.gz") %>%
  as.data.frame(.)
gene_df <- filter(ens_df, type == "gene", gene_biotype == "protein_coding")
gene_gr <- GRanges(seqnames = gene_df$seqnames,
                   ranges = IRanges(start = gene_df$start, end = gene_df$end),
                   gene_name = gene_df$gene_name)
names(gene_gr) <- gene_df$gene_name

# Read in the Manhattan input files
message("Reading in Manhattan input files...")
gwas_files <- list.files(input_dir, pattern = "manhattan-input_.*\\.txt$", full.names = TRUE)

if (length(gwas_files) == 0) {
  stop("No Manhattan input files found in the input directory: " %&% input_dir)
}

# Process each Manhattan input file
for (gwas_file in gwas_files) {
  gwas_name <- tools::file_path_sans_ext(basename(gwas_file))  # Extract the base name without extension
  
  message("Processing file: " %&% gwas_name)
  gwas.df <- fread(gwas_file)
  
  # No need to rename columns; use them directly as per your input file
  # Convert CHR to character and handle chromosome X
  gwas.df$CHR <- as.character(gwas.df$CHR)
  if ("23" %in% unique(gwas.df$CHR)){
    gwas.df$CHR <- gsub("23", gwas.df$CHR)
  }

  # Define the P value columns based on the structure of the current file
  p_columns <- c("P", "Pval_anc0", "Pval_anc1")

  # Process data frame
  plot.df <- process_gwas_df(gwas.df)

  # Generate Manhattan plots for each P column
  for (p_col in p_columns) {
    if (!(p_col %in% colnames(gwas.df))) {
      next  # Skip if the P column does not exist in the file
    }
    
    message("Generating Manhattan plot for: " %&% p_col %&% " in file: " %&% gwas_name)
    minp <- min(na.omit(gwas.df[[p_col]]))
    if (minp < 5e-8){
      plt.list <- manhat(plot.df, p_col, gene.gr = gene_gr)
    } else if (minp > 5e-8 & minp < 1e-6){
      plt.list <- manhat(plot.df, p_col, sig = 1e-6, gene.gr = gene_gr)
    } else {
      plt.list <- manhat(plot.df, p_col, sig = 1e-5, gene.gr = gene_gr)
    }
    
    m.plt <- plt.list[[1]]
    sig.df <- plt.list[[3]]
    
    # Save the plot and results
    ggsave(plot = m.plt, filename = plot_dir %&% gwas_name %&% "_" %&% p_col %&% ".png", width = 12, height = 6, type = "cairo-png")
    if (nrow(sig.df) > 0){
      write.table(x = sig.df, file = output_dir %&% "nearest-genes_" %&% gwas_name %&% "_" %&% p_col %&% ".txt",
                  quote = F, sep = "\t", col.names = T, row.names = F)
    }
  }
}

