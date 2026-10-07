# ============================================================
# Create patient metadata
# ============================================================

library(readxl)
library(dplyr)
library(stringr)

# 1. Input / output files

input_file <- file.path("metadata/original/INFO-PATIENTS_Methylseq.xlsx")

output_file <- file.path("metadata/patients_metadata.txt")

# 2. Read original patient information

patients <- read_excel(
  input_file,
  sheet = "Hoja 1",
  skip = 1
)

# 3. Build metadata table

patients_metadata <- patients %>%
  select(
    pathology = PATHOLOGY,
    gene_variant = `GENE VARIANT`,
    sample_id = `SAMPLE ID`,
    relationship = RELATIONSHIP,
    genetic_status = `GENETIC STATUS`,
    clinical_status = `CLINICAL STATUS`
  ) %>%
  mutate(
    patient_id = as.integer(str_remove(sample_id, "^P")),
    
    pathology = str_squish(pathology),
    gene_variant = str_squish(gene_variant),
    relationship = str_squish(relationship),
    genetic_status = str_squish(genetic_status),
    clinical_status = str_squish(clinical_status)
  ) %>%
  select(
    patient_id,
    pathology,
    gene_variant,
    relationship,
    genetic_status,
    clinical_status
  ) %>%
  mutate(
    across(
      everything(),
      ~ ifelse(is.na(.), "NA", .)
    )
  ) %>%
  arrange(patient_id)

# 4. Write tab-delimited TXT

write.table(
  patients_metadata,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

# 5. Confirmation

message("Metadata file created: ", output_file)
message("Number of patients: ", nrow(patients_metadata))