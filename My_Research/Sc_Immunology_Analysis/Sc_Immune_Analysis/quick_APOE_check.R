# Quick APOE Expression Check
# GSE222520 - Brain Tumor Leukocyte Single-Cell Analysis
# Efficient APOE analysis focusing on key patterns

library(dplyr)
library(ggplot2)

# Create output directories
dir.create("results/APOE_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/APOE_analysis", recursive = TRUE, showWarnings = FALSE)

# Load metadata
metadata <- read.csv("data/metadata/GSE222520_phenotype_data.csv", stringsAsFactors = FALSE)

# Create sample info
sample_info <- data.frame(
  sample_id = gsub(", replicate1, scRNAseq", "", metadata$title),
  geo_accession = metadata$geo_accession,
  tissue = metadata$`tissue:ch1`,
  diagnosis = metadata$`diagnosis:ch1`,
  genotype = metadata$`genotype:ch1`,
  primary_recurrent = metadata$`primary/recurrent:ch1`,
  treatment = metadata$`treatment:ch1`,
  stringsAsFactors = FALSE
)

# Define groups
sample_info$group <- case_when(
  grepl("^NGB", sample_info$sample_id) ~ "Normal",
  grepl("^IMP", sample_info$sample_id) ~ "Primary_IDH_Mutant",
  grepl("^IWP", sample_info$sample_id) ~ "Primary_IDH_WildType",
  grepl("^IMR", sample_info$sample_id) ~ "Recurrent_IDH_Mutant",
  grepl("^IWR", sample_info$sample_id) ~ "Recurrent_IDH_WildType"
)

sample_info$apoe_group <- case_when(
  sample_info$genotype == "IDH Mutant" ~ "IDH_Mutant",
  sample_info$genotype == "IDH Wild Type" ~ "IDH_WildType",
  TRUE ~ "Normal"
)

# Quick check of sample directories
sample_dirs <- list.dirs("data/raw/extracted", recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[!grepl("^\\.", sample_dirs)]

cat("Found", length(sample_dirs), "sample directories\n")

# Quick APOE gene check function
check_apoe_genes <- function(sample_dir) {
  sample_name <- gsub("GSM[0-9]+_", "", sample_dir)
  
  # Find matrix directory
  matrix_path <- file.path("data/raw/extracted", sample_dir)
  subdirs <- list.dirs(matrix_path, recursive = FALSE)
  if (length(subdirs) > 0) {
    matrix_path <- file.path(subdirs[1], "filtered_feature_bc_matrix")
  }
  
  if (!dir.exists(matrix_path)) {
    return(NULL)
  }
  
  # Read features file to check for APOE genes
  features_file <- file.path(matrix_path, "features.tsv.gz")
  if (!file.exists(features_file)) {
    return(NULL)
  }
  
  tryCatch({
    features <- readLines(gzfile(features_file))
    gene_names <- sapply(strsplit(features, "\t"), function(x) x[2])
    
    # Find APOE-related genes
    apoe_found <- gene_names[grepl("APO[EC]", gene_names, ignore.case = TRUE)]
    
    if (length(apoe_found) > 0) {
      return(data.frame(
        sample = sample_name,
        apoe_genes = paste(apoe_found, collapse = ", "),
        n_apoe_genes = length(apoe_found),
        stringsAsFactors = FALSE
      ))
    } else {
      return(data.frame(
        sample = sample_name,
        apoe_genes = "None",
        n_apoe_genes = 0,
        stringsAsFactors = FALSE
      ))
    }
  }, error = function(e) {
    return(data.frame(
      sample = sample_name,
      apoe_genes = "Error",
      n_apoe_genes = 0,
      stringsAsFactors = FALSE
    ))
  })
}

# Check all samples
cat("Checking APOE genes in all samples...\n")
apoe_results <- lapply(sample_dirs, check_apoe_genes)
apoe_results <- do.call(rbind, apoe_results)

# Combine with metadata
if (!is.null(apoe_results)) {
  apoe_results <- merge(apoe_results, sample_info, by.x = "sample", by.y = "sample_id", all.x = TRUE)
  
  # Save results
  write.csv(apoe_results, "results/APOE_analysis/APOE_gene_check.csv", row.names = FALSE)
  
  # Create summary
  cat("\n=== APOE GENE CHECK SUMMARY ===\n")
  cat("Total samples checked:", nrow(apoe_results), "\n")
  cat("Samples with APOE genes:", sum(apoe_results$n_apoe_genes > 0), "\n")
  
  # APOE genes found
  apoe_genes_found <- unique(unlist(strsplit(apoe_results$apoe_genes[apoe_results$apoe_genes != "None" & apoe_results$apoe_genes != "Error"], ", ")))
  cat("APOE-related genes found:", paste(apoe_genes_found, collapse = ", "), "\n")
  
  # Summary by group
  group_summary <- apoe_results %>%
    group_by(apoe_group) %>%
    summarise(
      total_samples = n(),
      samples_with_apoe = sum(n_apoe_genes > 0),
      apoe_rate = samples_with_apoe / total_samples * 100
    )
  
  print(group_summary)
  
  # Create visualization
  p1 <- ggplot(apoe_results, aes(x = apoe_group, fill = n_apoe_genes > 0)) +
    geom_bar() +
    theme_minimal() +
    labs(title = "APOE Gene Presence by IDH Genotype Group",
         x = "IDH Genotype Group", y = "Number of Samples", fill = "APOE Genes Present")
  
  ggsave("figures/APOE_analysis/APOE_gene_presence_by_group.pdf", p1, width = 10, height = 6)
  
  # Summary by diagnosis
  diagnosis_summary <- apoe_results %>%
    group_by(diagnosis) %>%
    summarise(
      total_samples = n(),
      samples_with_apoe = sum(n_apoe_genes > 0),
      apoe_rate = samples_with_apoe / total_samples * 100
    )
  
  print(diagnosis_summary)
  
  cat("\nResults saved to: results/APOE_analysis/APOE_gene_check.csv\n")
  cat("Figures saved to: figures/APOE_analysis/\n")
  
} else {
  cat("No results obtained from APOE gene check.\n")
}

cat("\n=== QUICK APOE CHECK COMPLETE ===\n")
