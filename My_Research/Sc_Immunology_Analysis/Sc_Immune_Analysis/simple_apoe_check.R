# Simple APOE Isoform Check
cat("=== SIMPLE APOE ISOFORM CHECK ===\n")

# Get sample directories
sample_dirs <- list.dirs("data/raw/extracted", recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[!grepl("^\\.", sample_dirs)]

cat("Found", length(sample_dirs), "sample directories\n")

# Check first sample
first_sample <- sample_dirs[1]
cat("Checking sample:", first_sample, "\n")

# Find matrix directory
matrix_path <- file.path("data/raw/extracted", first_sample)
subdirs <- list.dirs(matrix_path, recursive = FALSE)
if (length(subdirs) > 0) {
  matrix_path <- file.path(subdirs[1], "filtered_feature_bc_matrix")
}

if (dir.exists(matrix_path)) {
  # Read features file
  features_file <- file.path(matrix_path, "features.tsv.gz")
  if (file.exists(features_file)) {
    features <- readLines(gzfile(features_file))
    
    # Find APOE-related genes
    apoe_genes <- c()
    apoe_ensembl <- c()
    
    for (line in features) {
      parts <- strsplit(line, "\t")[[1]]
      gene_name <- parts[2]
      
      # Check for APOE-related genes
      if (grepl("APO[EC]", gene_name, ignore.case = TRUE)) {
        apoe_genes <- c(apoe_genes, gene_name)
        apoe_ensembl <- c(apoe_ensembl, parts[1])
      }
    }
    
    if (length(apoe_genes) > 0) {
      cat("\n=== APOE-RELATED GENES FOUND ===\n")
      cat("Total APOE-related genes:", length(apoe_genes), "\n\n")
      
      for (i in 1:length(apoe_genes)) {
        cat("Gene", i, ":\n")
        cat("  ENSEMBL ID:", apoe_ensembl[i], "\n")
        cat("  Gene Name:", apoe_genes[i], "\n")
        cat("\n")
      }
      
      # Check for specific APOE
      apoe_exact <- apoe_genes[apoe_genes == "APOE"]
      if (length(apoe_exact) > 0) {
        cat("=== APOE GENE FOUND ===\n")
        cat("APOE gene is present in the data\n")
        apoe_idx <- which(apoe_genes == "APOE")
        cat("ENSEMBL ID:", apoe_ensembl[apoe_idx], "\n")
      } else {
        cat("=== APOE GENE NOT FOUND ===\n")
        cat("APOE gene not found in exact match\n")
      }
      
      # Check for APOE variants
      apoe_variants <- apoe_genes[grepl("APOE", apoe_genes, ignore.case = TRUE)]
      if (length(apoe_variants) > 0) {
        cat("\n=== APOE VARIANTS FOUND ===\n")
        for (variant in apoe_variants) {
          cat("  ", variant, "\n")
        }
      }
      
      # Check chromosome 19 cluster
      chr19_genes <- c("APOE", "APOC1", "APOC2", "APOC4", "APOC4-APOC2")
      chr19_found <- apoe_genes[apoe_genes %in% chr19_genes]
      
      if (length(chr19_found) > 0) {
        cat("\n=== CHROMOSOME 19 APOE CLUSTER ===\n")
        cat("Chromosome 19 APOE cluster genes found:\n")
        for (gene in chr19_found) {
          cat("  ", gene, "\n")
        }
      }
      
      # Save results
      results <- data.frame(
        ensembl_id = apoe_ensembl,
        gene_name = apoe_genes,
        stringsAsFactors = FALSE
      )
      
      dir.create("results/APOE_analysis", recursive = TRUE, showWarnings = FALSE)
      write.csv(results, "results/APOE_analysis/apoe_genes_found.csv", row.names = FALSE)
      cat("\nAPOE genes list saved to: results/APOE_analysis/apoe_genes_found.csv\n")
      
    } else {
      cat("No APOE-related genes found in the data\n")
    }
  } else {
    cat("Features file not found\n")
  }
} else {
  cat("Matrix directory not found\n")
}

cat("\n=== APOE ISOFORM CHECK COMPLETE ===\n")
