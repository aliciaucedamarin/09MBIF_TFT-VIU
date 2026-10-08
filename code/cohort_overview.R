# ============================================================
# CHARACTERIZATION OF PATHOLOGY GROUPS
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Read metadata

metadata <- read_tsv(
  "metadata/patients_metadata.txt",
  na = c("NA"),
  show_col_types = FALSE
)

# 2. Create output directory

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. Variables represented within pathology

variables_to_plot <- c(
  "clinical_status",
  "gene_variant",
  "relationship",
  "genetic_status"
)

# 4. Figure titles

variable_labels <- c(
  clinical_status = "Clinical status",
  gene_variant = "Gene variant",
  relationship = "Relationship",
  genetic_status = "Genetic status"
)

# 5. Prepare data

pathology_characterization <- metadata %>%
  select(
    pathology,
    all_of(variables_to_plot)
  ) %>%
  mutate(
    pathology = ifelse(
      is.na(pathology),
      "NA",
      as.character(pathology)
    )
  ) %>%
  pivot_longer(
    cols = all_of(variables_to_plot),
    names_to = "variable",
    values_to = "category"
  ) %>%
  mutate(
    category = ifelse(
      is.na(category),
      "NA",
      as.character(category)
    )
  ) %>%
  count(
    pathology,
    variable,
    category,
    name = "n"
  ) %>%
  group_by(
    pathology,
    variable
  ) %>%
  mutate(
    percentage = 100 * n / sum(n)
  ) %>%
  ungroup()

# 6. Generate PDF

pdf(
  "results/figures/cohort_summary.pdf",
  width = 10,
  height = 7
)

# 7. Generate one figure per variable

for (var in variables_to_plot) {
  
  plot_data <- pathology_characterization %>%
    filter(variable == var)
  
  # Order pathology groups by total number of patients
  pathology_order <- plot_data %>%
    group_by(pathology) %>%
    summarise(
      total = sum(n),
      .groups = "drop"
    ) %>%
    arrange(desc(total)) %>%
    pull(pathology)
  
  plot_data <- plot_data %>%
    mutate(
      pathology = factor(
        pathology,
        levels = pathology_order
      )
    )
  
  # Labels
  
  plot_data <- plot_data %>%
    mutate(
      label = ifelse(
        percentage >= 0,
        paste0(
          round(percentage, 1),
          "%"
        ),
        ""
      )
    )
  
  # Plot
  
  p <- ggplot(
    plot_data,
    aes(
      x = pathology,
      y = percentage,
      fill = category
    )
  ) +
    
    # Stacked bars
    geom_col(
      width = 0.68,
      colour = "white",
      linewidth = 0.35
    ) +
    
    # Percentage labels
    geom_text(
      aes(
        label = label
      ),
      position = position_stack(
        vjust = 0.5
      ),
      size = 3.5,
      fontface = "bold",
      colour = "black"
    ) +
    
    # Y axis
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) {
        paste0(x, "%")
      },
      expand = expansion(
        mult = c(0, 0.01)
      )
    ) +
    
    # More contrasted palette
    scale_fill_manual(
      values = c(
        "#1B4F72",
        "#2874A6",
        "#5499C7",
        "#7FB3D5",
        "#148F77",
        "#45B39D",
        "#D68910",
        "#E67E22",
        "#A93226",
        "#7D3C98",
        "#566573"
      )
    ) +
    
    labs(
      title = variable_labels[var],
      subtitle = "Distribution within each pathology group",
      x = "Pathology",
      y = "Percentage of patients",
      fill = NULL
    ) +
    
    # Clean scientific style
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      
      # Title
      plot.title = element_text(
        size = 18,
        face = "bold",
        colour = "black",
        hjust = 0
      ),
      
      # Subtitle
      plot.subtitle = element_text(
        size = 11,
        colour = "grey35",
        hjust = 0,
        margin = margin(
          b = 15
        )
      ),
      
      # Axis titles
      axis.title.x = element_text(
        size = 11,
        face = "bold",
        colour = "black",
        margin = margin(
          t = 10
        )
      ),
      
      axis.title.y = element_text(
        size = 11,
        face = "bold",
        colour = "black",
        margin = margin(
          r = 10
        )
      ),
      
      # Axis labels
      axis.text.x = element_text(
        size = 10,
        colour = "black"
      ),
      
      axis.text.y = element_text(
        size = 10,
        colour = "black"
      ),
      
      # Legend
      legend.position = "right",
      
      legend.title = element_blank(),
      
      legend.text = element_text(
        size = 10,
        colour = "black"
      ),
      
      legend.key.height = unit(
        0.45,
        "cm"
      ),
      
      legend.key.width = unit(
        0.45,
        "cm"
      ),
      
      # IMPORTANT:
      # No background grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Clean panel
      panel.background = element_blank(),
      
      # White figure background
      plot.background = element_blank(),
      
      # Margins
      plot.margin = margin(
        15,
        20,
        15,
        15
      )
    )
  
  print(p)
}

dev.off()

# 8. Confirmation

cat("\n")
cat("============================================\n")
cat("PDF GENERATED\n")
cat("============================================\n")
cat("File:\n")
cat("results/figures/cohort_summary.pdf\n")
cat("============================================\n")