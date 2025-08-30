# Sample Data for Zijie's Great Mouse Management Software
# This file contains example data structures and functions

# Sample Cage Data (Enhanced Structure)
sample_cages <- data.frame(
  Cage_Name = c("Breeding-05 (03/07)", "Cx3-M1", "Tmem-M1", "Cx3-M2", "Tmem-M2", "Weaned-01"),
  Genotype = c("Cx3cr1(+/-)", "Cx3cr1(+/+)", "Tmem119(+/+)", "Cx3cr1(+/-)", "Tmem119(+/+)", "Cx3cr1(?/?)"),
  Sex = c("M", "M", "M", "F", "F", "Mixed"),
  Num_Mice = c(2, 1, 1, 2, 1, 4),
  Num_Pups = c(0, 0, 0, 0, 0, 4),
  Status = c("Breeding", "Holding", "Experiment", "Active", "Active", "Holding"),
  Parent_Cage = c("", "", "", "", "", "Breeding-05 (03/07)"),
  Genotyping_Done = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
  Genotype_Confirmed = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  Setup_Date = c("03/07", "01/15", "01/10", "01/20", "01/12", "01/16"),
  Notes = c("Breeding pair", "Control group", "Experiment group", "Young mice", "Health monitoring", "Weaned pups"),
  stringsAsFactors = FALSE
)

# Sample Mouse Data (Enhanced Structure)
sample_mice <- data.frame(
  Mouse_ID = c("Breeding-05-01", "Breeding-05-02", "Cx3-M1-01", "Tmem-M1-01", "Cx3-M2-01", "Cx3-M2-02", "W1-01", "W1-02", "W1-03", "W1-04"),
  Cage = c("Breeding-05 (03/07)", "Breeding-05 (03/07)", "Cx3-M1", "Tmem-M1", "Cx3-M2", "Cx3-M2", "Weaned-01", "Weaned-01", "Weaned-01", "Weaned-01"),
  Sex = c("M", "F", "M", "M", "F", "F", "M", "F", "M", "F"),
  DOB = c("2023-10-15", "2023-11-05", "2023-10-20", "2023-11-10", "2023-12-01", "2023-12-15", "2024-01-15", "2024-01-15", "2024-01-15", "2024-01-15"),
  Genotype = c("Cx3cr1(+/-)", "Cx3cr1(+/-)", "Cx3cr1(+/+)", "Tmem119(+/+)", "Cx3cr1(+/-)", "Cx3cr1(+/-)", "Cx3cr1(?/?)", "Cx3cr1(?/?)", "Cx3cr1(?/?)", "Cx3cr1(?/?)"),
  Status = c("Active", "Active", "Active", "Active", "Active", "Active", "Pup", "Pup", "Pup", "Pup"),
  Weight = c(25.5, 23.2, 26.8, 24.1, 22.5, 21.8, 18.2, 17.8, 19.1, 18.5),
  Health_Status = c("Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy"),
  Strain = c("C57BL/6J", "C57BL/6J", "C57BL/6J", "C57BL/6J", "C57BL/6J", "C57BL/6J", "C57BL/6J", "C57BL/6J", "C57BL/6J", "C57BL/6J"),
  Source = c("Jackson Labs", "Jackson Labs", "Internal", "Internal", "Internal", "Internal", "Internal", "Internal", "Internal", "Internal"),
  Lab_ID = c("LAB001", "LAB002", "LAB003", "LAB004", "LAB005", "LAB006", "LAB007", "LAB008", "LAB009", "LAB010"),
  Experiment_Group = c("Breeding", "Breeding", "Control", "Experiment", "Treatment", "Treatment", "Weaned", "Weaned", "Weaned", "Weaned"),
  Notes = c("Breeding male", "Breeding female", "Control group", "Experiment group", "Young female", "Young female", "Weaned", "Weaned", "Weaned", "Weaned"),
  Mother_ID = c(NA, NA, NA, NA, NA, NA, "Breeding-05-02", "Breeding-05-02", "Breeding-05-02", "Breeding-05-02"),
  Father_ID = c(NA, NA, NA, NA, NA, NA, "Breeding-05-01", "Breeding-05-01", "Breeding-05-01", "Breeding-05-01"),
  stringsAsFactors = FALSE
)

# Sample Gene Catalog Data
sample_genes <- data.frame(
  Gene = c("Cx3cr1", "Tmem119"),
  Symbols = c("+/+, +/-, -/-, WT, flox/+, flox/flox", "+/+, +/-, -/-, WT, KO, HET"),
  stringsAsFactors = FALSE
)

# Sample Breeding Data
sample_breeding <- data.frame(
  Breeding_ID = c("B001", "B002", "B003"),
  Male_ID = c("Breeding-05-01", "Cx3-M2-01", "Tmem-M1-01"),
  Female_ID = c("Breeding-05-02", "Cx3-M2-02", "Tmem-M1-02"),
  Start_Date = c("2024-01-01", "2024-01-05", "2024-01-10"),
  Expected_Birth = c("2024-01-21", "2024-01-25", "2024-01-30"),
  Status = c("Active", "Active", "Active"),
  Litter_Size = c(6, 5, 7),
  Notes = c("First breeding", "Second breeding", "Third breeding"),
  stringsAsFactors = FALSE
)

# Sample Health Records
sample_health <- data.frame(
  Record_ID = c("H001", "H002", "H003", "H004", "H005"),
  Mouse_ID = c("Breeding-05-01", "Breeding-05-02", "Cx3-M1-01", "Cx3-M2-01", "Cx3-M2-02"),
  Check_Date = c("2024-01-10", "2024-01-10", "2024-01-12", "2024-01-12", "2024-01-15"),
  Weight_g = c(25.5, 23.2, 26.8, 25.1, 24.3),
  Health_Status = c("Healthy", "Healthy", "Healthy", "Healthy", "Healthy"),
  Observations = c("Normal behavior", "Normal behavior", "Normal behavior", "Normal behavior", "Normal behavior"),
  Next_Check_Date = c("2024-01-24", "2024-01-24", "2024-01-26", "2024-01-26", "2024-01-29"),
  stringsAsFactors = FALSE
)

# Data Management Functions

# Function to get cage information
get_cage_info <- function(cage_name) {
  cage_data <- sample_cages[sample_cages$Cage_Name == cage_name, ]
  if (nrow(cage_data) == 0) {
    return(NULL)
  }
  return(cage_data)
}

# Function to get mice in a cage
get_mice_in_cage <- function(cage_name) {
  mice_data <- sample_mice[sample_mice$Cage == cage_name, ]
  return(mice_data)
}

# Function to get mouse information
get_mouse_info <- function(mouse_id) {
  mouse_data <- sample_mice[sample_mice$Mouse_ID == mouse_id, ]
  if (nrow(mouse_data) == 0) {
    return(NULL)
  }
  return(mouse_data)
}

# Function to get gene catalog
get_gene_catalog <- function() {
  return(sample_genes)
}

# Function to get breeding information
get_breeding_info <- function(mouse_id = NULL) {
  if (is.null(mouse_id)) {
    return(sample_breeding)
  } else {
    breeding_data <- sample_breeding[
      sample_breeding$Male_ID == mouse_id | sample_breeding$Female_ID == mouse_id, 
    ]
    return(breeding_data)
  }
}

# Function to get health records
get_health_records <- function(mouse_id = NULL) {
  if (is.null(mouse_id)) {
    return(sample_health)
  } else {
    health_data <- sample_health[sample_health$Mouse_ID == mouse_id, ]
    return(health_data)
  }
}

# Function to calculate cage statistics
calculate_cage_stats <- function(cage_name) {
  mice_in_cage <- get_mice_in_cage(cage_name)
  if (nrow(mice_in_cage) == 0) {
    return(list(
      total_mice = 0,
      males = 0,
      females = 0,
      pups = 0,
      avg_weight = 0
    ))
  }
  
  total_mice <- nrow(mice_in_cage)
  males <- sum(mice_in_cage$Sex == "M")
  females <- sum(mice_in_cage$Sex == "F")
  pups <- sum(mice_in_cage$Status == "Pup")
  avg_weight <- mean(mice_in_cage$Weight, na.rm = TRUE)
  
  return(list(
    total_mice = total_mice,
    males = males,
    females = females,
    pups = pups,
    avg_weight = avg_weight
  ))
}

# Function to get parent information
get_parent_info <- function(mouse_id) {
  mouse_data <- get_mouse_info(mouse_id)
  if (is.null(mouse_data)) {
    return(NULL)
  }
  
  mother_id <- mouse_data$Mother_ID[1]
  father_id <- mouse_data$Father_ID[1]
  
  mother_info <- if (!is.na(mother_id)) get_mouse_info(mother_id) else NULL
  father_info <- if (!is.na(father_id)) get_mouse_info(father_id) else NULL
  
  return(list(
    mother = mother_info,
    father = father_info
  ))
}

# Function to get offspring
get_offspring <- function(mouse_id) {
  offspring <- sample_mice[
    sample_mice$Mother_ID == mouse_id | sample_mice$Father_ID == mouse_id, 
  ]
  return(offspring)
}

# Function to search mice by criteria
search_mice <- function(criteria = list()) {
  results <- sample_mice
  
  if (!is.null(criteria$cage)) {
    results <- results[results$Cage == criteria$cage, ]
  }
  
  if (!is.null(criteria$sex)) {
    results <- results[results$Sex == criteria$sex, ]
  }
  
  if (!is.null(criteria$status)) {
    results <- results[results$Status == criteria$status, ]
  }
  
  if (!is.null(criteria$genotype)) {
    results <- results[grepl(criteria$genotype, results$Genotype), ]
  }
  
  if (!is.null(criteria$search_text)) {
    search_pattern <- tolower(criteria$search_text)
    results <- results[
      grepl(search_pattern, tolower(results$Mouse_ID)) |
      grepl(search_pattern, tolower(results$Notes)) |
      grepl(search_pattern, tolower(results$Genotype)), 
    ]
  }
  
  return(results)
}

# Function to search cages by criteria
search_cages <- function(criteria = list()) {
  results <- sample_cages
  
  if (!is.null(criteria$status)) {
    results <- results[results$Status == criteria$status, ]
  }
  
  if (!is.null(criteria$genotyping_done)) {
    results <- results[results$Genotyping_Done == criteria$genotyping_done, ]
  }
  
  if (!is.null(criteria$has_pups)) {
    if (criteria$has_pups) {
      results <- results[results$Num_Pups > 0, ]
    } else {
      results <- results[results$Num_Pups == 0, ]
    }
  }
  
  if (!is.null(criteria$search_text)) {
    search_pattern <- tolower(criteria$search_text)
    results <- results[grepl(search_pattern, tolower(results$Cage_Name)), ]
  }
  
  return(results)
}

# Export functions for data
export_cage_data <- function(cage_name, format = "csv") {
  cage_info <- get_cage_info(cage_name)
  mice_in_cage <- get_mice_in_cage(cage_name)
  
  if (format == "csv") {
    # Return as list for CSV export
    return(list(
      cage_info = cage_info,
      mice = mice_in_cage
    ))
  } else {
    # Return as data frame for other formats
    return(list(
      cage_info = cage_info,
      mice = mice_in_cage
    ))
  }
}

export_all_data <- function(format = "csv") {
  if (format == "csv") {
    return(list(
      cages = sample_cages,
      mice = sample_mice,
      genes = sample_genes,
      breeding = sample_breeding,
      health = sample_health
    ))
  } else {
    return(list(
      cages = sample_cages,
      mice = sample_mice,
      genes = sample_genes,
      breeding = sample_breeding,
      health = sample_health
    ))
  }
}

# Print summary of sample data
cat("Sample Data Summary for Zijie's Great Mouse Management Software:\n")
cat("• Cages:", nrow(sample_cages), "\n")
cat("• Mice:", nrow(sample_mice), "\n")
cat("• Genes:", nrow(sample_genes), "\n")
cat("• Breeding Records:", nrow(sample_breeding), "\n")
cat("• Health Records:", nrow(sample_health), "\n")
cat("\nEnhanced Features Available:\n")
cat("• Parent tracking (Mother_ID, Father_ID)\n")
cat("• Gene catalog with genotype symbols\n")
cat("• Comprehensive mouse data (weight, health, strain, etc.)\n")
cat("• Accurate pup counting and status tracking\n")
cat("• Search and filter functions\n")
cat("• Export capabilities\n")
