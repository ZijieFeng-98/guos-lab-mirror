# Targeted APOE Analysis
# Using exact genes found in the data

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
  tissue = metadata$tissue.ch1,
  diagnosis = metadata$diagnosis.ch1,
  genotype = metadata$genotype.ch1,
  primary_recurrent = metadata$primary.recurrent.ch1,
  stringsAsFactors = FALSE
)

# Define groups
sample_info$apoe_group <- case_when(
  sample_info$genotype == "IDH Mutant" ~ "IDH_Mutant",
  sample_info$genotype == "IDH Wild Type" ~ "IDH_WildType",
  TRUE ~ "Normal"
)

# Exact APOE genes found in the data
apoe_genes <- c("APOE", "APOC1", "APOC2", "APOC3", "APOC4", "APOC4-APOC2")
apoe_ensembl <- c("ENSG00000130203", "ENSG00000130208", "ENSG00000234906", 
                  "ENSG00000110245", "ENSG00000267467", "ENSG00000224916")

# Function to extract APOE expression
extract_apoe_expression <- function(sample_dir) {
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
    
    # Find APOE gene indices
    apoe_indices <- c()
    apoe_names <- c()
    
    for (i in 1:length(features)) {
      parts <- strsplit(features[i], "\t")[[1]]
      if (parts[2] %in% apoe_genes) {
        apoe_indices <- c(apoe_indices, i)
        apoe_names <- c(apoe_names, parts[2])
      }
    }
    
    if (length(apoe_indices) > 0) {
      # Read matrix file
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
      
      # Calculate expression for each APOE gene
      apoe_results <- data.frame()
      for (i in 1:length(apoe_indices)) {
        gene_idx <- apoe_indices[i]
        gene_name <- apoe_names[i]
        
        # Get expression for this gene
        gene_entries <- parsed_entries[parsed_entries[, 1] == gene_idx, , drop = FALSE]
        
        if (nrow(gene_entries) > 0) {
          mean_expr <- mean(gene_entries[, 3])
          n_expressing_cells <- nrow(gene_entries)
          pct_expressing <- n_expressing_cells / n_cells * 100
          max_expr <- max(gene_entries[, 3])
          median_expr <- median(gene_entries[, 3])
        } else {
          mean_expr <- 0
          n_expressing_cells <- 0
          pct_expressing <- 0
          max_expr <- 0
          median_expr <- 0
        }
        
        apoe_results <- rbind(apoe_results, data.frame(
          sample = sample_name,
          gene = gene_name,
          mean_expression = mean_expr,
          median_expression = median_expr,
          max_expression = max_expr,
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
        median_expression = 0,
        max_expression = 0,
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
      median_expression = 0,
      max_expression = 0,
      n_expressing_cells = 0,
      total_cells = 0,
      pct_expressing = 0,
      stringsAsFactors = FALSE
    ))
  })
}

# Get sample directories
sample_dirs <- list.dirs("data/raw/extracted", recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[!grepl("^\\.", sample_dirs)]

cat("=== TARGETED APOE ANALYSIS ===\n")
cat("Targeting genes:", paste(apoe_genes, collapse = ", "), "\n")
cat("Analyzing", length(sample_dirs), "samples...\n")

# Extract APOE data from all samples
all_results <- lapply(sample_dirs, extract_apoe_expression)
all_results <- do.call(rbind, all_results)

# Combine with metadata
if (!is.null(all_results) && nrow(all_results) > 0) {
  all_results <- merge(all_results, sample_info, by.x = "sample", by.y = "sample_id", all.x = TRUE)
  
  # Save results
  write.csv(all_results, "results/APOE_analysis/targeted_APOE_results.csv", row.names = FALSE)
  
  # Create summary
  cat("\n=== TARGETED APOE SUMMARY ===\n")
  cat("Total samples analyzed:", length(unique(all_results$sample)), "\n")
  
  # APOE expression summary by group
  apoe_summary <- all_results %>%
    filter(gene == "APOE") %>%
    group_by(apoe_group) %>%
    summarise(
      n_samples = n(),
      mean_expression = mean(mean_expression, na.rm = TRUE),
      median_expression = median(median_expression, na.rm = TRUE),
      mean_pct_expressing = mean(pct_expressing, na.rm = TRUE),
      max_expression = max(max_expression, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\nAPOE Expression by IDH Genotype Group:\n")
  print(apoe_summary)
  
  # All APOE family summary
  family_summary <- all_results %>%
    filter(gene %in% apoe_genes) %>%
    group_by(gene, apoe_group) %>%
    summarise(
      n_samples = n(),
      mean_expression = mean(mean_expression, na.rm = TRUE),
      mean_pct_expressing = mean(pct_expressing, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\nAPOE Family Expression Summary:\n")
  print(family_summary)
  
  # Create visualizations
  if (nrow(all_results) > 0) {
    # APOE expression by group
    p1 <- ggplot(all_results %>% filter(gene == "APOE"), 
                 aes(x = apoe_group, y = mean_expression, fill = apoe_group)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      labs(title = "APOE Expression by IDH Genotype Group",
           x = "IDH Genotype Group", y = "Mean Expression", fill = "Group")
    
    ggsave("figures/APOE_analysis/targeted_APOE_expression.pdf", p1, width = 10, height = 6)
    
    # All APOE family genes
    p2 <- ggplot(all_results %>% filter(gene %in% apoe_genes), 
                 aes(x = gene, y = mean_expression, fill = apoe_group)) +
      geom_boxplot() +
      theme_minimal() +
      labs(title = "APOE Family Gene Expression",
           x = "Gene", y = "Mean Expression", fill = "Group")
    
    ggsave("figures/APOE_analysis/targeted_APOE_family.pdf", p2, width = 12, height = 6)
    
    # Expression percentage
    p3 <- ggplot(all_results %>% filter(gene == "APOE"), 
                 aes(x = apoe_group, y = pct_expressing, fill = apoe_group)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      labs(title = "APOE Expression Percentage",
           x = "IDH Genotype Group", y = "Percentage of Expressing Cells", fill = "Group")
    
    ggsave("figures/APOE_analysis/targeted_APOE_percentage.pdf", p3, width = 10, height = 6)
  }
  
  cat("\nResults saved to: results/APOE_analysis/targeted_APOE_results.csv\n")
  cat("Figures saved to: figures/APOE_analysis/\n")
  
} else {
  cat("No results obtained from APOE analysis.\n")
}

cat("\n=== TARGETED APOE ANALYSIS COMPLETE ===\n")
