# ============================================================
# COMPLETE MATRIX ANNOTATION
# ============================================================

# 1. LOAD LIBRARIES
# Make sure the hg19 TxDb package is installed.
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# install.packages(c("readxl", "openxlsx", "tidyverse"))

library(readxl)
library(openxlsx)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(tidyverse)

# 2. LOAD THE ORIGINAL DATA FROM EXCEL

input_file <- file.path(
  "data/original/all_24_content_clean_annot_nodup.bed.xlsx"
)

raw_data <- read_excel(
  input_file,
  sheet = "all_24_content_clean_annot_nodu",
  skip = 0
)

# Preserve a copy of the original dataset for integrity checks
raw_data_original <- raw_data

# 3. CONVERT THE COMPLETE DATASET INTO A GRANGES OBJECT
# Process all rows to match their genomic coordinates against hg19.

all_regions_gr <- GRanges(
  seqnames = raw_data$`#CHROM`,
  ranges = IRanges(
    start = raw_data$START,
    end = raw_data$END
  )
)

# 4. PERFORM GLOBAL GENOMIC ANNOTATION USING HG19

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

annotated_peaks_all <- annotatePeak(
  all_regions_gr,
  tssRegion = c(-2000, 2000),
  TxDb = txdb_hg19,
  annoDb = "org.Hs.eg.db"
)

# Convert annotation results into a structured data frame

annotation_temp_df <- as.data.frame(annotated_peaks_all)

global_annotation_df <- data.frame(
  `#CHROM` = as.character(annotation_temp_df$seqnames),
  START = annotation_temp_df$start,
  END = annotation_temp_df$end,
  NEW_GENE = annotation_temp_df$SYMBOL,
  REGIONAL_TYPE = annotation_temp_df$annotation,
  DISTANCE_TO_TSS = annotation_temp_df$distanceToTSS,
  check.names = FALSE
)

# 5. IDENTIFY MULTIPLE ANNOTATIONS PER GENOMIC COORDINATE

# Count the number of annotation records for each coordinate combination.
# Multiple records may correspond to different genes or transcripts.

duplicated_annotations <- global_annotation_df %>%
  dplyr::count(
    `#CHROM`, START, END,
    name = "N_ANNOTATIONS"
  ) %>%
  dplyr::filter(N_ANNOTATIONS > 1)

cat(
  "Coordinate combinations with multiple annotations:",
  nrow(duplicated_annotations),
  "\n"
)

# 6. CONSOLIDATE ANNOTATIONS TO PREVENT ROW MULTIPLICATION

# Keep one row per coordinate combination.
# Preserve distinct gene symbols, regional annotations and TSS distances
# by concatenating multiple values with semicolons.

global_annotation_unique <- global_annotation_df %>%
  dplyr::group_by(`#CHROM`, START, END) %>%
  dplyr::summarise(
    NEW_GENE = {
      x <- as.character(NEW_GENE)
      x <- unique(x[!is.na(x) & x != "" & x != "."])
      paste(x, collapse = ";")
    },
    REGIONAL_TYPE = {
      x <- as.character(REGIONAL_TYPE)
      x <- unique(x[!is.na(x) & x != ""])
      paste(x, collapse = ";")
    },
    DISTANCE_TO_TSS = {
      x <- DISTANCE_TO_TSS
      x <- unique(x[!is.na(x)])
      paste(x, collapse = ";")
    },
    .groups = "drop"
  )

# Verify that the annotation table now has unique coordinate combinations

stopifnot(
  !anyDuplicated(
    global_annotation_unique[c("#CHROM", "START", "END")]
  )
)

# 7. INTEGRATE ANNOTATIONS WHILE PRESERVING THE ORIGINAL ROWS

# Add the annotations using a many-to-one join.
# Fill missing or placeholder gene names with mapped hg19 gene symbols.
# Preserve valid gene names already present in the original dataset.

raw_data_annotated_final <- raw_data %>%
  dplyr::left_join(
    global_annotation_unique,
    by = c("#CHROM", "START", "END"),
    relationship = "many-to-one"
  ) %>%
  dplyr::mutate(
    GEN = dplyr::if_else(
      is.na(GEN) | GEN %in% c(".", "-1"),
      dplyr::na_if(NEW_GENE, ""),
      as.character(GEN)
    )
  ) %>%
  dplyr::select(-NEW_GENE) %>%
  dplyr::relocate(
    REGIONAL_TYPE,
    DISTANCE_TO_TSS,
    .after = GEN
  )

# 8. VERIFY DATASET DIMENSIONS AND ROW PRESERVATION

cat("\n===== DATASET DIMENSIONS =====\n")

cat("Original rows:", nrow(raw_data_original), "\n")
cat("Final rows:", nrow(raw_data_annotated_final), "\n")
cat("Original columns:", ncol(raw_data_original), "\n")
cat("Final columns:", ncol(raw_data_annotated_final), "\n")

# Stop if the number of rows has unexpectedly changed

stopifnot(
  nrow(raw_data_annotated_final) == nrow(raw_data_original)
)

cat("PASS: The original number of rows has been preserved.\n")

# 9. VERIFY THAT ORIGINAL COLUMNS REMAIN UNCHANGED, EXCEPT FOR GEN

# Check all original columns except GEN, which is intentionally updated.

original_columns <- setdiff(
  names(raw_data_original),
  "GEN"
)

missing_columns <- setdiff(
  original_columns,
  names(raw_data_annotated_final)
)

if (length(missing_columns) > 0) {
  cat("FAIL: Original columns missing from the final dataset:\n")
  print(missing_columns)
} else {
  
  unchanged_columns <- vapply(
    original_columns,
    function(col) {
      isTRUE(
        all.equal(
          raw_data_original[[col]],
          raw_data_annotated_final[[col]],
          check.attributes = FALSE
        )
      )
    },
    logical(1)
  )
  
  if (all(unchanged_columns)) {
    cat(
      "PASS: All original columns except GEN remain unchanged.\n"
    )
  } else {
    cat("WARNING: Some original columns have changed:\n")
    print(names(unchanged_columns)[!unchanged_columns])
  }
}

# 10. VERIFY THAT VALID ORIGINAL GENE NAMES WERE PRESERVED

# Compare gene names for rows that had valid original values.
# This comparison is valid because the join preserves row order
# and the row count has been verified.

original_genes <- as.character(raw_data_original$GEN)
final_genes <- as.character(raw_data_annotated_final$GEN)

valid_original_genes <- !is.na(original_genes) &
  !original_genes %in% c(".", "-1", "")

genes_preserved <- all(
  original_genes[valid_original_genes] ==
    final_genes[valid_original_genes]
)

cat(
  "Original valid gene names preserved:",
  genes_preserved,
  "\n"
)

# 11. REPORT ANNOTATION COVERAGE

cat("\n===== ANNOTATION SUMMARY =====\n")

cat(
  "Rows with a regional annotation:",
  sum(!is.na(raw_data_annotated_final$REGIONAL_TYPE)),
  "\n"
)

cat(
  "Rows with a TSS distance:",
  sum(!is.na(raw_data_annotated_final$DISTANCE_TO_TSS)),
  "\n"
)

# 12. SAVE AND EXPORT THE ANNOTATED DATASET TO EXCEL

output_file <- "data/raw_data_annotated_hg19.xlsx"

write.xlsx(
  raw_data_annotated_final,
  file = output_file,
  sheetName = "all_24_content_annotated",
  overwrite = TRUE
)

cat("\nAnnotated dataset exported to:", output_file, "\n")

# 13. FINAL INTEGRITY CHECK

# Confirm that the final dataset has the expected dimensions
# and that the original data columns were preserved.

stopifnot(
  nrow(raw_data_annotated_final) == nrow(raw_data_original),
  all(original_columns %in% names(raw_data_annotated_final)),
  all(unchanged_columns),
  genes_preserved
)

cat(
  "\nPASS: All integrity checks completed successfully.\n"
)