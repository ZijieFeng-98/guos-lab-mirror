# Fast APOE Analysis
# GSE222520 - Brain Tumor Leukocyte Single-Cell Analysis

library(dplyr)
library(ggplot2)

# Create output directories
dir.create("results/APOE_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/APOE_analysis", recursive = TRUE, showWarnings = FALSE)

# Load metadata
metadata <- read.csv("data/metadata/GSE222520_phenotype_data.csv", stringsAsFactors = FALSE)

# Create sample info with correct column names
sample_info <- data.frame(
  sample_id = gsub(", replicate1, scRNAseq", "", metadata$title),
  geo_accession = metadata$geo_accession,
  tissue = metadata$tissue.ch1,
  diagnosis = metadata$diagnosis.ch1,
  genotype = metadata$genotype.ch1,
  primary_recurrent = metadata$primary.recurrent.ch1,
  treatment = metadata$treatment.ch1,
  cell_type = metadata$cell.type.ch1,
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

# Get sample directories
sample_dirs <- list.dirs("data/raw/extracted", recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[!grepl("^\\.", sample_dirs)]

# APOE genes to check
apoe_genes <- c("APOE", "APOC1", "APOC2", "APOC3", "APOC4")

# Function to check APOE expression in one sample
check_apoe_expression <- function(sample_dir) {
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
  
  # Read features file
  features_file <- file.path(matrix_path, "features.tsv.gz")
  if (!file.exists(features_file)) {
    return(NULL)
  }
  
  tryCatch({
    features <- readLines(gzfile(features_file))
    gene_names <- sapply(strsplit(features, "\t"), function(x) x[2])
    
    # Find APOE genes
    apoe_indices <- which(gene_names %in% apoe_genes)
    
    if (length(apoe_indices) > 0) {
      # Read matrix file to get expression
      matrix_file <- file.path(matrix_path, "matrix.mtx.gz")
      matrix_lines <- readLines(gzfile(matrix_file))
      
      # Parse header
      header <- strsplit(matrix_lines[1], " ")[[1]]
      n_genes <- as.numeric(header[1])
      n_cells <- as.numeric(header[2])
      
      # Parse matrix entries
      entries <- matrix_lines[2:length(matrix_lines)]
      parsed_entries <- do.call(rbind, lapply(entries, function(x) {
        parts <- strsplit(x, " ")[[1]]
        c(as.numeric(parts[1]), as.numeric(parts[2]), as.numeric(parts[3]))
      }))
      
      # Calculate APOE expression
      apoe_results <- data.frame()
      for (i in seq_along(apoe_indices)) {
        gene_idx <- apoe_indices[i]
        gene_name <- gene_names[gene_idx]
        
        # Get expression for this gene
        gene_entries <- parsed_entries[parsed_entries[, 1] == gene_idx, , drop = FALSE]
        
        if (nrow(gene_entries) > 0) {
          mean_expr <- mean(gene_entries[, 3])
          n_expressing_cells <- nrow(gene_entries)
          pct_expressing <- n_expressing_cells / n_cells * 100
        } else {
          mean_expr <- 0
          n_expressing_cells <- 0
          pct_expressing <- 0
        }
        
        apoe_results <- rbind(apoe_results, data.frame(
          sample = sample_name,
          gene = gene_name,
          mean_expression = mean_expr,
          n_expressing_cells = n_expressing_cells,
          total_cells = n_cells,
          pct_expressing = pct_expressing,
          stringsAsFactors = FALSE
        ))
      }
      
      return(apoe_results)
    } else {
      return(data.frame(
        sample = sample_name,
        gene = "None",
        mean_expression = 0,
        n_expressing_cells = 0,
        total_cells = 0,
        pct_expressing = 0,
        stringsAsFactors = FALSE
      ))
    }
  }, error = function(e) {
    return(data.frame(
      sample = sample_name,
      gene = "Error",
      mean_expression = 0,
      n_expressing_cells = 0,
      total_cells = 0,
      pct_expressing = 0,
      stringsAsFactors = FALSE
    ))
  })
}

# Check all samples
cat("Checking APOE expression in all samples...\n")
all_results <- lapply(sample_dirs, check_apoe_expression)
all_results <- do.call(rbind, all_results)

# Combine with metadata
if (!is.null(all_results) && nrow(all_results) > 0) {
  all_results <- merge(all_results, sample_info, by.x = "sample", by.y = "sample_id", all.x = TRUE)
  
  # Save results
  write.csv(all_results, "results/APOE_analysis/APOE_expression_results.csv", row.names = FALSE)
  
  # Create summary
  cat("\n=== APOE EXPRESSION SUMMARY ===\n")
  cat("Total samples analyzed:", length(unique(all_results$sample)), "\n")
  
  # APOE expression summary
  apoe_summary <- all_results %>%
    filter(gene == "APOE") %>%
    group_by(apoe_group) %>%
    summarise(
      n_samples = n(),
      mean_expression = mean(mean_expression, na.rm = TRUE),
      mean_pct_expressing = mean(pct_expressing, na.rm = TRUE),
      .groups = 'drop'
    )
  
  print(apoe_summary)
  
  # Create visualizations
  if (nrow(all_results) > 0) {
    # APOE expression by group
    p1 <- ggplot(all_results %>% filter(gene == "APOE"), 
                 aes(x = apoe_group, y = mean_expression, fill = apoe_group)) +
      geom_boxplot() +
      theme_minimal() +
      labs(title = "APOE Expression by IDH Genotype Group",
           x = "IDH Genotype Group", y = "Mean Expression", fill = "Group")
    
    ggsave("figures/APOE_analysis/APOE_expression_by_group.pdf", p1, width = 10, height = 6)
    
    # APOE percentage expressing cells
    p2 <- ggplot(all_results %>% filter(gene == "APOE"), 
                 aes(x = apoe_group, y = pct_expressing, fill = apoe_group)) +
      geom_boxplot() +
      theme_minimal() +
      labs(title = "APOE Expression Percentage by IDH Genotype Group",
           x = "IDH Genotype Group", y = "Percentage of Expressing Cells", fill = "Group")
    
    ggsave("figures/APOE_analysis/APOE_pct_expressing_by_group.pdf", p2, width = 10, height = 6)
    
    # All APOE family genes
    p3 <- ggplot(all_results %>% filter(gene %in% apoe_genes), 
                 aes(x = gene, y = mean_expression, fill = apoe_group)) +
      geom_boxplot() +
      theme_minimal() +
      labs(title = "APOE Family Gene Expression by IDH Genotype Group",
           x = "Gene", y = "Mean Expression", fill = "Group")
    
    ggsave("figures/APOE_analysis/APOE_family_expression.pdf", p3, width = 12, height = 6)
  }
  
  cat("\nResults saved to: results/APOE_analysis/APOE_expression_results.csv\n")
  cat("Figures saved to: figures/APOE_analysis/\n")
  
} else {
  cat("No results obtained from APOE expression analysis.\n")
}

cat("\n=== FAST APOE ANALYSIS COMPLETE ===\n")
