# ==============================================================================
# PIPELINE: PCA -> LIMMA -> VOLCANO PLOT -> CANDIDATE DMR EXTRACTION
# ==============================================================================

# 1. INSTALL AND LOAD REQUIRED PACKAGES
if (!require("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!require("readxl", quietly = TRUE)) install.packages("readxl")
if (!require("ggfortify", quietly = TRUE)) install.packages("ggfortify")
if (!require("limma", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("limma")
}
if (!require("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")

library(tidyverse)
library(readxl)
library(ggfortify)
library(limma)
library(EnhancedVolcano)


# ==============================================================================
# 2. DATA IMPORTATION
# ==============================================================================

input_file <- file.path(
  "data/raw_data_annotated_hg19.xlsx"
)

raw_data <- read_excel(
  input_file,
  sheet = "all_24_content_annotated",
  skip = 0
) 

metadata <- read_tsv(
  "metadata/patients_metadata.txt",
  na = c("NA"),
  show_col_types = FALSE
)


# ==============================================================================
# 3. DATA PREPARATION AND SAMPLE ALIGNMENT
# ==============================================================================

metadata <- metadata %>%
  mutate(
    Muestra_Col = paste0(patient_id, "_Norm")
  ) %>%
  filter(
    Muestra_Col %in% colnames(bed_matrix)
  )

norm_matrix <- bed_matrix %>%
  select(all_of(metadata$Muestra_Col))

norm_matrix_clean <- norm_matrix %>%
  mutate(
    across(
      everything(),
      ~ as.numeric(as.character(.))
    )
  )

numeric_matrix <- as.matrix(norm_matrix_clean)

unique_regions <- make.unique(
  bed_matrix$`REGION ID`
)

rownames(numeric_matrix) <- unique_regions

annotation <- bed_matrix %>%
  mutate(
    Unique_ID = unique_regions
  ) %>%
  select(
    Unique_ID,
    `#CHROM`,
    START,
    END,
    `REGION ID`,
    chr,
    gene_start,
    gene_end,
    GEN
  )

group <- factor(
  metadata$clinical_status,
  levels = c("HEALTHY", "AFFECTED")
)


# ==============================================================================
# 4. DATA TRANSFORMATION
# ==============================================================================

# Log2 transformation with pseudocount
log_matrix <- log2(numeric_matrix + 1)


# ==============================================================================
# 5. PRINCIPAL COMPONENT ANALYSIS
# ==============================================================================

row_variances <- apply(
  log_matrix,
  1,
  var
)

log_matrix_filtered <- log_matrix[
  row_variances > 0,
]

pca_results <- prcomp(
  t(log_matrix_filtered),
  scale. = TRUE
)

variance_explained <- summary(
  pca_results
)$importance[2, ] * 100

pc1_var <- round(
  variance_explained[1],
  2
)

pc2_var <- round(
  variance_explained[2],
  2
)

pca_plot <- autoplot(
  pca_results,
  data = metadata,
  colour = "clinical_status",
  size = 4
) +
  theme_minimal() +
  labs(
    title = "PCA: Global Sample Clustering",
    x = paste0(
      "Principal Component 1 (",
      pc1_var,
      "%)"
    ),
    y = paste0(
      "Principal Component 2 (",
      pc2_var,
      "%)"
    ),
    colour = "Clinical Status"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

print(pca_plot)


# ==============================================================================
# 6. DIFFERENTIAL ANALYSIS WITH LIMMA
# ==============================================================================

design <- model.matrix(
  ~ group
)

fit <- lmFit(
  log_matrix,
  design
)

fit <- eBayes(
  fit
)

stats_table <- topTable(
  fit,
  coef = "groupAFFECTED",
  number = Inf,
  adjust.method = "BH"
)


# ==============================================================================
# 7. VOLCANO PLOT
# ==============================================================================

volcano_plot <- EnhancedVolcano(
  stats_table,
  lab = rownames(stats_table),
  x = "logFC",
  y = "P.Value",
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.5,
  labSize = 3.0,
  title = "Volcano Plot: AFFECTED vs HEALTHY",
  subtitle = "Differential analysis using limma",
  legendPosition = "bottom"
)

print(volcano_plot)


# ==============================================================================
# 8. CANDIDATE DMR EXTRACTION
# ==============================================================================

final_results <- stats_table %>%
  rownames_to_column(
    var = "Unique_ID"
  ) %>%
  left_join(
    annotation,
    by = "Unique_ID"
  )

candidate_DMRs <- final_results %>%
  filter(
    P.Value < 0.05,
    abs(logFC) >= 1.0
  ) %>%
  arrange(
    P.Value
  )

write_csv(
  candidate_DMRs,
  "results/tables/Significant_DMRs_Report.csv"
)

print("====================================================================")
print(
  paste(
    "Process complete! Found",
    nrow(candidate_DMRs),
    "candidate DMRs (P < 0.05, |logFC| >= 1.0)."
  )
)
print("====================================================================")