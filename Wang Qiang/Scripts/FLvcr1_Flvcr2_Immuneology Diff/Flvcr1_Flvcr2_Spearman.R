# ============================
# FLVCR1/FLVCR2 correlations by cohort → Genes as rows (readable labels)
# ============================

# install.packages(c("readxl","dplyr","stringr","pheatmap","RColorBrewer")) # if needed
library(readxl)
library(dplyr)
library(stringr)
library(pheatmap)
library(RColorBrewer)

path <- "~/Desktop/Guo's lab/Wang Qiang/Raw_Data/Flvcr1 and Flvcr2.xlsx"   # adjust to absolute path if needed
raw  <- read_excel(path, sheet = 1, col_names = FALSE)

to_num <- function(x) suppressWarnings(as.numeric(x))

# ---- 1) Find cohort title rows (contain TCGA/GDC labels) ----
is_cohort_row <- apply(raw, 1, function(r) any(stringr::str_detect(as.character(r), "(TCGA|GDC)"), na.rm = TRUE))
cohort_idx <- which(is_cohort_row)

get_title <- function(i){
  row <- as.character(unlist(raw[i, ]))
  row <- row[!is.na(row) & trimws(row) != ""]
  if (length(row) == 0) "Unknown" else row[1]
}
titles <- sapply(cohort_idx, get_title)

blocks <- data.frame(
  start = cohort_idx,
  end   = c(cohort_idx[-1] - 1, nrow(raw)),
  title = titles,
  stringsAsFactors = FALSE
)

# ---- 2) Keep the 3 cohorts of interest ----
want_patterns <- c("GBM", "Pancreatic", "Ovarian")  # edit to swap cohorts
blocks$keep <- sapply(blocks$title, function(tt) any(str_detect(tt, str_c(want_patterns, collapse = "|"))))
blocks_sel <- blocks %>% filter(keep)
if (nrow(blocks_sel) == 0) stop("No matching cohorts found. Check 'want_patterns' or Excel contents.")

# ---- 3) Parse one cohort block into (gene, FLVCR1, FLVCR2) table ----
parse_block <- function(raw, start_idx, end_idx) {
  hdr_row <- NA_integer_
  for (i in seq(start_idx, end_idx)) {
    rowv <- as.character(unlist(raw[i, ]))
    if (any(rowv == "FLVCR1", na.rm = TRUE) && any(rowv == "FLVCR2", na.rm = TRUE)) {
      hdr_row <- i; break
    }
  }
  if (is.na(hdr_row)) return(NULL)
  
  row_hdr <- as.character(unlist(raw[hdr_row, ]))
  spearman_cols <- which(row_hdr == "Spearman's correlation")
  if (length(spearman_cols) < 2) return(NULL)
  
  gene_cols <- spearman_cols - 1
  data_rows <- seq(hdr_row + 1, end_idx)
  
  left  <- raw[data_rows, c(gene_cols[1], spearman_cols[1])]
  right <- raw[data_rows, c(gene_cols[2], spearman_cols[2])]
  names(left)  <- c("gene",  "FLVCR1")
  names(right) <- c("gene2", "FLVCR2")
  
  left <- left %>%
    mutate(gene = trimws(as.character(gene)),
           FLVCR1 = to_num(FLVCR1)) %>%
    filter(!is.na(gene), gene != "", !is.na(FLVCR1))
  
  right <- right %>%
    mutate(gene2 = trimws(as.character(gene2)),
           FLVCR2 = to_num(FLVCR2)) %>%
    filter(!is.na(gene2), gene2 != "", !is.na(FLVCR2))
  
  # collapse duplicates within the cohort (mean)
  left_c  <- left  %>% group_by(gene)  %>% summarise(FLVCR1 = mean(FLVCR1, na.rm = TRUE), .groups = "drop")
  right_c <- right %>% group_by(gene2) %>% summarise(FLVCR2 = mean(FLVCR2, na.rm = TRUE), .groups = "drop")
  
  dat <- full_join(left_c, right_c, by = c("gene" = "gene2")) %>%
    filter(!(is.na(FLVCR1) & is.na(FLVCR2))) %>%
    distinct(gene, .keep_all = TRUE)
  
  dat
}

# ---- 4) Parse selected cohorts ----
cohort_tables <- list()
cohort_labels <- c()
for (k in seq_len(nrow(blocks_sel))) {
  b <- blocks_sel[k, ]
  tab <- parse_block(raw, b$start, b$end)
  if (!is.null(tab) && nrow(tab) > 0) {
    cohort_tables[[length(cohort_tables) + 1]] <- tab
    cohort_labels <- c(cohort_labels, b$title)
  }
}
if (length(cohort_tables) == 0) stop("Matched cohorts exist but parsing failed.")

# ---- 5) Within-block column ordering (cluster each block's genes) ----
order_block <- function(tab) {
  m <- rbind(FLVCR1 = tab$FLVCR1, FLVCR2 = tab$FLVCR2)
  colnames(m) <- tab$gene
  d  <- dist(t(m))
  hc <- hclust(d, method = "average")
  ord <- hc$order
  list(mat = m[, ord, drop = FALSE], ord = ord)
}
block_mats <- lapply(cohort_tables, order_block)

# ---- 6) Make colnames unique by appending cohort tags; store clean labels ----
tags <- sapply(cohort_labels, function(t) {
  if (grepl("GBM", t, ignore.case = TRUE)) "GBM"
  else if (grepl("Pancre", t, ignore.case = TRUE)) "PAAD"
  else if (grepl("Ovar", t, ignore.case = TRUE)) "OV"
  else gsub("[^A-Za-z0-9]+", "", substr(t, 1, 8))
})
display_labels_list <- vector("list", length(block_mats))
for (i in seq_along(block_mats)) {
  genes_i <- colnames(block_mats[[i]]$mat)
  display_labels_list[[i]] <- genes_i
  colnames(block_mats[[i]]$mat) <- make.unique(paste0(genes_i, "|", tags[i]))
}
display_labels <- unlist(display_labels_list, use.names = FALSE)

# ---- 7) Combine blocks + gaps + annotations (for original orientation) ----
mat_combined <- do.call(cbind, lapply(block_mats, function(x) x$mat))
block_sizes  <- sapply(block_mats, function(x) ncol(x$mat))
gaps         <- cumsum(block_sizes); gaps <- gaps[-length(gaps)]

CancerType <- rep(cohort_labels, times = block_sizes)
sign_flip  <- sign(mat_combined["FLVCR1", ]) != sign(mat_combined["FLVCR2", ])

# ============================
# ---- 8) Flip so GENES ARE ROWS (readable labels) ----
# ============================
mat_t <- t(mat_combined)  # rows = genes (with cohort-tagged IDs), cols = FLVCR1/FLVCR2

# Row annotations (one per gene)
row_annot <- data.frame(
  CancerType = CancerType,
  SignFlip   = ifelse(sign_flip, "Opposite", "Same"),
  row.names  = rownames(mat_t)
)

# Column annotation (only 2 columns)
col_annot <- data.frame(Variable = c("FLVCR1", "FLVCR2"))
rownames(col_annot) <- colnames(mat_t)

# Colors
ct_levels  <- unique(CancerType)
ct_palette <- setNames(brewer.pal(max(3, length(ct_levels)), "Set2")[seq_along(ct_levels)], ct_levels)
ann_colors_rows <- list(
  CancerType = ct_palette,
  SignFlip   = c(Opposite = "#d95f02", Same = "#1b9e77")
)
ann_colors_cols <- list(
  Variable = c(FLVCR1 = "#4575b4", FLVCR2 = "#d73027")
)

# Row gaps to separate the 3 cohorts
gaps_row <- cumsum(block_sizes); if (length(gaps_row) > 1) gaps_row <- gaps_row[-length(gaps_row)] else gaps_row <- NULL

# Show clean gene names (no cohort tag)
labels_row <- display_labels

# Color map
pal    <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(101)
breaks <- seq(-1, 1, length.out = 101)

# ---- 9) Plot tall heatmap (genes as rows) ----
# PNG
png("FLVCR1_FLVCR2_GenesAsRows.png", width = 1800, height = 3000, res = 250)
pheatmap(
  mat_t,
  cluster_rows = FALSE,          # keep within-block order we computed
  cluster_cols = FALSE,          # only two columns
  gaps_row = gaps_row,           # 3 cohort blocks
  color = pal, breaks = breaks,
  border_color = NA,
  main = "Spearman correlations: Genes (rows) vs FLVCR1/FLVCR2 (cols)",
  annotation_row = row_annot,
  annotation_col = col_annot,
  annotation_colors = c(ann_colors_rows, ann_colors_cols),
  show_rownames = TRUE,
  labels_row = labels_row,
  fontsize_row = 9,              # increase/decrease for readability
  cellheight = 14,               # increase for more spacing per gene
  fontsize_col = 14
)
dev.off()

# PDF (vector)
pdf("FLVCR1_FLVCR2_GenesAsRows.pdf", width = 8, height = 14)
pheatmap(
  mat_t,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gaps_row,
  color = pal, breaks = breaks,
  border_color = NA,
  main = "Spearman correlations: Genes (rows) vs FLVCR1/FLVCR2 (cols)",
  annotation_row = row_annot,
  annotation_col = col_annot,
  annotation_colors = c(ann_colors_rows, ann_colors_cols),
  show_rownames = TRUE,
  labels_row = labels_row,
  fontsize_row = 8.5,
  cellheight = 12,
  fontsize_col = 12
)
dev.off()

cat("\nCohorts:", paste(unique(CancerType), collapse=" | "),
    "\nGenes total:", nrow(mat_t),
    "\nSign flips:", sum(row_annot$SignFlip == "Opposite", na.rm = TRUE), "/", nrow(mat_t), "\n")
