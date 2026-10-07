# ============================================================
# COHORT DESCRIPTION - FIGURES
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Read the metadata file

metadata <- read_tsv(
  "metadata/patients_metadata.txt",
  na = c("NA"),
  show_col_types = FALSE
)

# 2. Create the output directory

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. General graphical style 

theme_tfm <- function() {
  
  theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        size = 17,
        face = "bold",
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 11,
        color = "grey30",
        hjust = 0
      ),
      axis.title = element_text(
        size = 12,
        face = "bold"
      ),
      axis.text = element_text(
        size = 11,
        color = "black"
      ),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5
      ),
      legend.title = element_text(
        size = 11,
        face = "bold"
      ),
      legend.text = element_text(
        size = 10
      ),
      legend.position = "right",
      plot.margin = margin(
        15, 20, 15, 15
      )
    )
}

# 4. Clinical status

clinical_summary <- metadata %>%
  count(clinical_status) %>%
  mutate(
    percentage = 100 * n / sum(n)
  )


p_clinical <- ggplot(
  clinical_summary,
  aes(
    x = clinical_status,
    y = n
  )
) +
  geom_col(
    fill = "steelblue",
    width = 0.65
  ) +
  geom_text(
    aes(
      label = paste0(
        n,
        " (",
        round(percentage, 1),
        "%)"
      )
    ),
    vjust = -0.5,
    size = 4.2
  ) +
  labs(
    title = "Clinical status of the cohort",
    x = "Clinical status",
    y = "Number of patients"
  ) +
  expand_limits(
    y = max(clinical_summary$n) * 1.15
  ) +
  theme_tfm()


# 5. Pathology

pathology_summary <- metadata %>%
  count(pathology) %>%
  mutate(
    percentage = 100 * n / sum(n)
  )


p_pathology <- ggplot(
  pathology_summary,
  aes(
    x = reorder(pathology, n),
    y = n
  )
) +
  geom_col(
    fill = "steelblue",
    width = 0.65
  ) +
  geom_text(
    aes(
      label = paste0(
        n,
        " (",
        round(percentage, 1),
        "%)"
      )
    ),
    hjust = -0.1,
    size = 4
  ) +
  coord_flip() +
  labs(
    title = "Pathological diagnosis",
    x = "Pathology",
    y = "Number of patients"
  ) +
  expand_limits(
    y = max(pathology_summary$n) * 1.15
  ) +
  theme_tfm()

# 6. Gene variant

gene_summary <- metadata %>%
  count(gene_variant) %>%
  mutate(
    percentage = 100 * n / sum(n)
  )


p_gene <- ggplot(
  gene_summary,
  aes(
    x = reorder(gene_variant, n),
    y = n
  )
) +
  geom_col(
    fill = "steelblue",
    width = 0.65
  ) +
  geom_text(
    aes(
      label = paste0(
        n,
        " (",
        round(percentage, 1),
        "%)"
      )
    ),
    hjust = -0.1,
    size = 4
  ) +
  coord_flip() +
  labs(
    title = "Gene variant distribution",
    x = "Gene variant",
    y = "Number of patients"
  ) +
  expand_limits(
    y = max(gene_summary$n) * 1.15
  ) +
  theme_tfm()

# 7. Genetic status

genetic_summary <- metadata %>%
  count(genetic_status) %>%
  mutate(
    percentage = 100 * n / sum(n)
  )


p_genetic <- ggplot(
  genetic_summary,
  aes(
    x = genetic_status,
    y = n
  )
) +
  geom_col(
    fill = "steelblue",
    width = 0.65
  ) +
  geom_text(
    aes(
      label = paste0(
        n,
        " (",
        round(percentage, 1),
        "%)"
      )
    ),
    vjust = -0.5,
    size = 4.2
  ) +
  labs(
    title = "Genetic status of the cohort",
    x = "Genetic status",
    y = "Number of patients"
  ) +
  expand_limits(
    y = max(genetic_summary$n) * 1.15
  ) +
  theme_tfm()

# 8. Relationship

relationship_summary <- metadata %>%
  count(relationship) %>%
  mutate(
    percentage = 100 * n / sum(n)
  )


p_relationship <- ggplot(
  relationship_summary,
  aes(
    x = reorder(relationship, n),
    y = n
  )
) +
  geom_col(
    fill = "steelblue",
    width = 0.65
  ) +
  geom_text(
    aes(
      label = paste0(
        n,
        " (",
        round(percentage, 1),
        "%)"
      )
    ),
    hjust = -0.1,
    size = 4
  ) +
  coord_flip() +
  labs(
    title = "Relationship distribution",
    x = "Relationship",
    y = "Number of patients"
  ) +
  expand_limits(
    y = max(relationship_summary$n) * 1.15
  ) +
  theme_tfm()

# 9. Pathology × Clinical status

pathology_clinical <- metadata %>%
  count(
    pathology,
    clinical_status
  ) %>%
  group_by(pathology) %>%
  mutate(
    percentage_within_pathology =
      100 * n / sum(n)
  ) %>%
  ungroup()


p_pathology_clinical <- ggplot(
  pathology_clinical,
  aes(
    x = pathology,
    y = percentage_within_pathology,
    fill = clinical_status
  )
) +
  geom_col(
    width = 0.7
  ) +
  labs(
    title = "Clinical status by pathology",
    x = "Pathology",
    y = "Percentage within pathology",
    fill = "Clinical status"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  theme_tfm() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# 10. Gene variant × Clinical status

gene_clinical <- metadata %>%
  count(
    gene_variant,
    clinical_status
  ) %>%
  group_by(gene_variant) %>%
  mutate(
    percentage_within_gene =
      100 * n / sum(n)
  ) %>%
  ungroup()


p_gene_clinical <- ggplot(
  gene_clinical,
  aes(
    x = gene_variant,
    y = percentage_within_gene,
    fill = clinical_status
  )
) +
  geom_col(
    width = 0.7
  ) +
  labs(
    title = "Clinical status by gene variant",
    x = "Gene variant",
    y = "Percentage within gene variant",
    fill = "Clinical status"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  theme_tfm() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# 11. Genetic status × Clinical status

genetic_clinical <- metadata %>%
  count(
    genetic_status,
    clinical_status
  ) %>%
  group_by(genetic_status) %>%
  mutate(
    percentage_within_genetic_status =
      100 * n / sum(n)
  ) %>%
  ungroup()


p_genetic_clinical <- ggplot(
  genetic_clinical,
  aes(
    x = genetic_status,
    y = percentage_within_genetic_status,
    fill = clinical_status
  )
) +
  geom_col(
    width = 0.7
  ) +
  labs(
    title = "Clinical status by genetic status",
    x = "Genetic status",
    y = "Percentage within genetic status",
    fill = "Clinical status"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  theme_tfm()

# 12. Missing data summary

missing_summary <- metadata %>%
  summarise(
    across(
      everything(),
      ~ sum(is.na(.))
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "missing_n"
  ) %>%
  mutate(
    total_n = nrow(metadata),
    missing_percentage =
      100 * missing_n / total_n
  )


p_missing <- ggplot(
  missing_summary,
  aes(
    x = reorder(
      variable,
      missing_percentage
    ),
    y = missing_percentage
  )
) +
  geom_col(
    fill = "steelblue",
    width = 0.65
  ) +
  geom_text(
    aes(
      label = paste0(
        round(
          missing_percentage,
          1
        ),
        "%"
      )
    ),
    hjust = -0.1,
    size = 3.8
  ) +
  coord_flip() +
  labs(
    title = "Missing data across cohort variables",
    x = "Variable",
    y = "Missing data"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%")
  ) +
  expand_limits(
    y = max(
      missing_summary$missing_percentage
    ) * 1.15
  ) +
  theme_tfm()

# 13. Clinical group sizes and imbalance ratio

group_sizes <- metadata %>%
  count(clinical_status) %>%
  mutate(
    percentage =
      100 * n / sum(n),
    imbalance_ratio =
      max(n) / n
  )


p_group_sizes <- ggplot(
  group_sizes,
  aes(
    x = clinical_status,
    y = n
  )
) +
  geom_col(
    fill = "steelblue",
    width = 0.65
  ) +
  geom_text(
    aes(
      label = paste0(
        n,
        " (",
        round(percentage, 1),
        "%)"
      )
    ),
    vjust = -0.5,
    size = 4.2
  ) +
  labs(
    title = "Clinical group sizes",
    x = "Clinical status",
    y = "Number of patients"
  ) +
  expand_limits(
    y = max(group_sizes$n) * 1.15
  ) +
  theme_tfm()

# 14. Save all figures to a single PDF file

pdf(
  "results/figures/cohort_description.pdf",
  width = 10,
  height = 7
)


# Figure 1
print(p_clinical)

# Figure 2
print(p_pathology)

# Figure 3
print(p_gene)

# Figure 4
print(p_genetic)

# Figure 5
print(p_relationship)

# Figure 6
print(p_pathology_clinical)

# Figure 7
print(p_gene_clinical)

# Figure 8
print(p_genetic_clinical)

# Figure 9
print(p_missing)

# Figure 10
print(p_group_sizes)


# Close PDF
dev.off()

# 15. Save summary tables for categorical variables

categorical_vars <- c(
  "pathology",
  "gene_variant",
  "relationship",
  "genetic_status",
  "clinical_status"
)


describe_categorical <- function(
    data,
    variable
) {
  
  data %>%
    count(
      .data[[variable]],
      name = "n"
    ) %>%
    mutate(
      variable = variable,
      category = ifelse(
        is.na(.data[[variable]]),
        "NA",
        as.character(.data[[variable]])
      ),
      percentage =
        round(
          100 * n / sum(n),
          1
        )
    ) %>%
    select(
      variable = .data[["variable"]],
      category = .data[["category"]],
      n = .data[["n"]],
      percentage = .data[["percentage"]]
    )
}


cohort_summary <- bind_rows(
  lapply(
    categorical_vars,
    function(x) {
      describe_categorical(
        metadata,
        x
      )
    }
  )
)


write_tsv(
  cohort_summary,
  "results/tables/cohort_categorical_summary.txt"
)


write_tsv(
  missing_summary,
  "results/tables/cohort_missing_values.txt"
)


write_tsv(
  group_sizes,
  "results/tables/clinical_group_sizes.txt"
)

# 16. Final messages

cat("\n")
cat("============================================\n")
cat("COHORT FIGURES COMPLETED\n")
cat("============================================\n")
cat("Patients:", nrow(metadata), "\n")
cat("Variables:", ncol(metadata), "\n")
cat("\n")
cat("PDF generated:\n")
cat("results/figures/cohort_description.pdf\n")
cat("============================================\n")