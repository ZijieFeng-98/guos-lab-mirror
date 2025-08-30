# Check Cell Types in Single-Cell Dataset
# Determine if dataset contains normal cells vs tumor cells

library(dplyr)
library(ggplot2)

cat("=== CHECKING CELL TYPES IN SINGLE-CELL DATASET ===\n")

# Load metadata
metadata <- read.csv("data/metadata/GSE222520_phenotype_data.csv", stringsAsFactors = FALSE)

# Create comprehensive sample info
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

cat("Total samples in dataset:", nrow(sample_info), "\n\n")

# Check tissue types
cat("=== TISSUE TYPES ===\n")
tissue_summary <- sample_info %>%
  group_by(tissue) %>%
  summarise(
    n_samples = n(),
    .groups = 'drop'
  )
print(tissue_summary)

# Check diagnoses
cat("\n=== DIAGNOSES ===\n")
diagnosis_summary <- sample_info %>%
  group_by(diagnosis) %>%
  summarise(
    n_samples = n(),
    .groups = 'drop'
  )
print(diagnosis_summary)

# Check cell types
cat("\n=== CELL TYPES ===\n")
celltype_summary <- sample_info %>%
  group_by(cell_type) %>%
  summarise(
    n_samples = n(),
    .groups = 'drop'
  )
print(celltype_summary)

# Define normal vs tumor samples
sample_info$sample_category <- case_when(
  grepl("^NGB", sample_info$sample_id) ~ "Normal_Brain",
  grepl("^IMP|^IWP|^IMR|^IWR", sample_info$sample_id) ~ "Tumor"
)

sample_info$detailed_category <- case_when(
  grepl("^NGB", sample_info$sample_id) ~ "Normal_Brain",
  grepl("^IMP", sample_info$sample_id) ~ "Primary_IDH_Mutant",
  grepl("^IWP", sample_info$sample_id) ~ "Primary_IDH_WildType",
  grepl("^IMR", sample_info$sample_id) ~ "Recurrent_IDH_Mutant",
  grepl("^IWR", sample_info$sample_id) ~ "Recurrent_IDH_WildType"
)

# Summary by sample category
cat("\n=== SAMPLE CATEGORIES ===\n")
category_summary <- sample_info %>%
  group_by(sample_category) %>%
  summarise(
    n_samples = n(),
    .groups = 'drop'
  )
print(category_summary)

cat("\n=== DETAILED SAMPLE CATEGORIES ===\n")
detailed_summary <- sample_info %>%
  group_by(detailed_category) %>%
  summarise(
    n_samples = n(),
    .groups = 'drop'
  )
print(detailed_summary)

# Check for normal samples specifically
normal_samples <- sample_info[sample_info$sample_category == "Normal_Brain", ]
tumor_samples <- sample_info[sample_info$sample_category == "Tumor", ]

cat("\n=== NORMAL SAMPLES DETAILS ===\n")
if (nrow(normal_samples) > 0) {
  cat("Found", nrow(normal_samples), "normal brain samples:\n")
  for (i in 1:nrow(normal_samples)) {
    cat("  ", normal_samples$sample_id[i], ":", normal_samples$diagnosis[i], "\n")
  }
} else {
  cat("No normal brain samples found!\n")
}

cat("\n=== TUMOR SAMPLES DETAILS ===\n")
cat("Found", nrow(tumor_samples), "tumor samples:\n")
tumor_diagnoses <- tumor_samples %>%
  group_by(diagnosis) %>%
  summarise(
    n_samples = n(),
    .groups = 'drop'
  )
print(tumor_diagnoses)

# Check sample naming patterns
cat("\n=== SAMPLE NAMING PATTERNS ===\n")
sample_patterns <- sample_info %>%
  mutate(
    pattern = case_when(
      grepl("^NGB", sample_id) ~ "NGB (Normal Brain)",
      grepl("^IMP", sample_id) ~ "IMP (Primary IDH Mutant)",
      grepl("^IWP", sample_id) ~ "IWP (Primary IDH Wild Type)",
      grepl("^IMR", sample_id) ~ "IMR (Recurrent IDH Mutant)",
      grepl("^IWR", sample_id) ~ "IWR (Recurrent IDH Wild Type)",
      TRUE ~ "Unknown"
    )
  ) %>%
  group_by(pattern) %>%
  summarise(
    n_samples = n(),
    .groups = 'drop'
  )
print(sample_patterns)

# Create visualizations
if (nrow(sample_info) > 0) {
  # Sample distribution by category
  p1 <- ggplot(sample_info, aes(x = sample_category, fill = sample_category)) +
    geom_bar() +
    theme_minimal() +
    labs(title = "Sample Distribution: Normal vs Tumor",
         x = "Sample Category", y = "Number of Samples", fill = "Category")
  
  ggsave("figures/sample_distribution.pdf", p1, width = 8, height = 6)
  
  # Detailed sample distribution
  p2 <- ggplot(sample_info, aes(x = detailed_category, fill = detailed_category)) +
    geom_bar() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Detailed Sample Distribution",
         x = "Sample Category", y = "Number of Samples", fill = "Category")
  
  ggsave("figures/detailed_sample_distribution.pdf", p2, width = 10, height = 6)
  
  # Diagnosis distribution
  p3 <- ggplot(sample_info, aes(x = diagnosis, fill = sample_category)) +
    geom_bar() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Diagnosis Distribution by Sample Category",
         x = "Diagnosis", y = "Number of Samples", fill = "Sample Category")
  
  ggsave("figures/diagnosis_distribution.pdf", p3, width = 12, height = 6)
}

# Save summary results
write.csv(sample_info, "results/sample_cell_type_summary.csv", row.names = FALSE)

# Final summary
cat("\n=== FINAL SUMMARY ===\n")
if (nrow(normal_samples) > 0) {
  cat("✅ Dataset CONTAINS normal cells\n")
  cat("   - Normal brain samples:", nrow(normal_samples), "\n")
  cat("   - Tumor samples:", nrow(tumor_samples), "\n")
  cat("   - Normal samples:", paste(normal_samples$sample_id, collapse = ", "), "\n")
} else {
  cat("❌ Dataset contains ONLY tumor cells\n")
  cat("   - No normal brain samples found\n")
  cat("   - All samples are from tumor tissue\n")
}

cat("\n=== RECOMMENDATIONS ===\n")
if (nrow(normal_samples) > 0) {
  cat("✅ You can compare normal vs tumor APOE expression\n")
  cat("✅ Use normal samples as controls for tumor analysis\n")
  cat("✅ Consider separate analyses for normal and tumor groups\n")
} else {
  cat("⚠️  No normal controls available\n")
  cat("⚠️  Can only analyze APOE expression within tumor samples\n")
  cat("⚠️  Consider comparing different tumor types/stages instead\n")
}

cat("\nResults saved to: results/sample_cell_type_summary.csv\n")
cat("Figures saved to: figures/\n")

cat("\n=== CELL TYPE CHECK COMPLETE ===\n")
