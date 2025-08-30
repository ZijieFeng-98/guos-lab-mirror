# Simple data structure check
cat("=== SIMPLE DATA CHECK ===\n")

# Check metadata
cat("1. Checking metadata...\n")
metadata <- read.csv("data/metadata/GSE222520_phenotype_data.csv", stringsAsFactors = FALSE)
cat("Metadata columns:", paste(colnames(metadata), collapse = ", "), "\n")
cat("Number of samples:", nrow(metadata), "\n")

# Check sample directories
cat("\n2. Checking sample directories...\n")
sample_dirs <- list.dirs("data/raw/extracted", recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[!grepl("^\\.", sample_dirs)]
cat("Sample directories found:", length(sample_dirs), "\n")
cat("First few:", paste(head(sample_dirs, 3), collapse = ", "), "\n")

# Check one sample structure
if (length(sample_dirs) > 0) {
  cat("\n3. Checking structure of first sample...\n")
  first_sample <- sample_dirs[1]
  cat("First sample:", first_sample, "\n")
  
  # List contents
  sample_path <- file.path("data/raw/extracted", first_sample)
  subdirs <- list.dirs(sample_path, recursive = FALSE)
  cat("Subdirectories:", paste(subdirs, collapse = ", "), "\n")
  
  if (length(subdirs) > 0) {
    matrix_path <- file.path(subdirs[1], "filtered_feature_bc_matrix")
    cat("Matrix path:", matrix_path, "\n")
    cat("Matrix exists:", dir.exists(matrix_path), "\n")
    
    if (dir.exists(matrix_path)) {
      files <- list.files(matrix_path)
      cat("Matrix files:", paste(files, collapse = ", "), "\n")
      
      # Try to read features file
      features_file <- file.path(matrix_path, "features.tsv.gz")
      if (file.exists(features_file)) {
        cat("Features file exists\n")
        # Read first few lines
        features <- readLines(gzfile(features_file), n = 5)
        cat("First 5 features:", paste(features, collapse = "\n"), "\n")
        
        # Check for APOE
        all_features <- readLines(gzfile(features_file))
        gene_names <- sapply(strsplit(all_features, "\t"), function(x) x[2])
        apoe_genes <- gene_names[grepl("APO[EC]", gene_names, ignore.case = TRUE)]
        cat("APOE-related genes found:", paste(apoe_genes, collapse = ", "), "\n")
      }
    }
  }
}

cat("\n=== CHECK COMPLETE ===\n")
