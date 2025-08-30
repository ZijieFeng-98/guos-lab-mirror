# Check APOE Isoforms and Gene Names
# First identify what's actually in the data

cat("=== CHECKING APOE ISOFORMS AND GENE NAMES ===\n")

# Get sample directories
sample_dirs <- list.dirs("data/raw/extracted", recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[!grepl("^\\.", sample_dirs)]

cat("Found", length(sample_dirs), "sample directories\n")

# Check first sample to identify APOE-related genes
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
    
    # Parse gene information
    gene_info <- do.call(rbind, lapply(features, function(x) {
      parts <- strsplit(x, "\t")[[1]]
      c(parts[1], parts[2], parts[3])  # ENSEMBL_ID, Gene_Name, Type
    }))
    
    colnames(gene_info) <- c("ensembl_id", "gene_name", "type")
    
    # Find all APOE-related genes
    apoe_patterns <- c("APOE", "APOC", "APOA", "APOB", "APOD", "APOF", "APOH", "APOL")
    apoe_genes <- c()
    
    for (pattern in apoe_patterns) {
      matches <- gene_info[grepl(pattern, gene_info[, "gene_name"], ignore.case = TRUE), ]
      if (nrow(matches) > 0) {
        apoe_genes <- rbind(apoe_genes, matches)
      }
    }
    
    if (nrow(apoe_genes) > 0) {
      cat("\n=== APOE-RELATED GENES FOUND ===\n")
      cat("Total APOE-related genes:", nrow(apoe_genes), "\n\n")
      
      # Display each gene with details
      for (i in 1:nrow(apoe_genes)) {
        cat("Gene", i, ":\n")
        cat("  ENSEMBL ID:", apoe_genes[i, "ensembl_id"], "\n")
        cat("  Gene Name:", apoe_genes[i, "gene_name"], "\n")
        cat("  Type:", apoe_genes[i, "type"], "\n")
        cat("\n")
      }
      
      # Check for specific APOE isoforms
      cat("=== APOE ISOFORM CHECK ===\n")
      apoe_exact <- gene_info[gene_info[, "gene_name"] == "APOE", ]
      if (nrow(apoe_exact) > 0) {
        cat("APOE gene found:\n")
        cat("  ENSEMBL ID:", apoe_exact[1, "ensembl_id"], "\n")
        cat("  Gene Name:", apoe_exact[1, "gene_name"], "\n")
        cat("  Type:", apoe_exact[1, "type"], "\n")
      } else {
        cat("APOE gene not found in exact match\n")
      }
      
      # Check for APOE variants/isoforms
      apoe_variants <- gene_info[grepl("APOE", gene_info[, "gene_name"], ignore.case = TRUE), ]
      if (nrow(apoe_variants) > 0) {
        cat("\nAPOE variants found:\n")
        for (i in 1:nrow(apoe_variants)) {
          cat("  ", apoe_variants[i, "gene_name"], " (", apoe_variants[i, "ensembl_id"], ")\n")
        }
      }
      
      # Check chromosome 19 genes
      cat("\n=== CHROMOSOME 19 APOE CLUSTER ===\n")
      # APOE gene cluster on chr19: 44,905,791-44,910,393
      chr19_apoe_genes <- c("APOE", "APOC1", "APOC2", "APOC4", "APOC4-APOC2")
      chr19_found <- gene_info[gene_info[, "gene_name"] %in% chr19_apoe_genes, ]
      
      if (nrow(chr19_found) > 0) {
        cat("Chromosome 19 APOE cluster genes found:\n")
        for (i in 1:nrow(chr19_found)) {
          cat("  ", chr19_found[i, "gene_name"], " (", chr19_found[i, "ensembl_id"], ")\n")
        }
      } else {
        cat("No chromosome 19 APOE cluster genes found\n")
      }
      
      # Save gene list for analysis
      write.csv(apoe_genes, "results/APOE_analysis/apoe_genes_found.csv", row.names = FALSE)
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
