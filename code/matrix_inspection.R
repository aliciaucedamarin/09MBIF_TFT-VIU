
# ============================================================
# ANALISIS EXPLORATORIO USANDO SOLO COLUMNAS _Norm
# Grupos: HEALTHY vs AFFECTED
# ============================================================


# 1. INSTALAR Y CARGAR PAQUETES --------------------------------

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "limma",
  "readr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(
  lapply(required_packages, library, character.only = TRUE)
)


# 2. RUTAS Y CONFIGURACION --------------------------------------

excel_file <- "data/raw_data_annotated_hg19.xlsx"
metadata_file <- "metadata/patients_metadata.txt"
excel_sheet <- "all_24_content_annotated"

output_tables <- "results/tables"
output_figures <- "results/figures"

dir.create(output_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figures, recursive = TRUE, showWarnings = FALSE)

fdr_threshold <- 0.05
top_regions_heatmap <- 50


# 3. COMPROBAR ARCHIVOS -----------------------------------------

if (!file.exists(excel_file)) {
  stop("No se encuentra el archivo Excel: ", excel_file)
}

if (!file.exists(metadata_file)) {
  stop("No se encuentra el archivo de metadatos: ", metadata_file)
}


# 4. LEER DATOS -------------------------------------------------

message("Leyendo archivo Excel...")

raw_data <- readxl::read_excel(
  path = excel_file,
  sheet = excel_sheet
) |>
  as.data.frame(check.names = FALSE)

message("Leyendo metadatos...")

metadata <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A")
)

# Limpiar espacios en los nombres de las columnas
names(raw_data) <- trimws(names(raw_data))
names(metadata) <- trimws(names(metadata))

# Columnas que deben existir en los metadatos
metadata_columns <- c(
  "patient_id",
  "pathology",
  "gene_variant",
  "relationship",
  "genetic_status",
  "clinical_status"
)

missing_metadata_columns <- setdiff(
  metadata_columns,
  names(metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "Faltan estas columnas en los metadatos: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

metadata$patient_id <- trimws(as.character(metadata$patient_id))
metadata$clinical_status <- toupper(
  trimws(as.character(metadata$clinical_status))
)

if (any(is.na(metadata$patient_id) | metadata$patient_id == "")) {
  stop("Hay patient_id vacíos en los metadatos.")
}

if (anyDuplicated(metadata$patient_id)) {
  duplicated_ids <- unique(
    metadata$patient_id[duplicated(metadata$patient_id)]
  )
  
  stop(
    "Hay patient_id duplicados en los metadatos: ",
    paste(duplicated_ids, collapse = ", ")
  )
}


# 5. IDENTIFICAR EXCLUSIVAMENTE LAS COLUMNAS _Norm -------------

norm_columns <- grep(
  "_Norm$",
  names(raw_data),
  value = TRUE
)

if (length(norm_columns) < 2) {
  stop(
    "Se necesitan al menos dos columnas terminadas en '_Norm'. ",
    "Comprueba los nombres de las columnas del Excel."
  )
}

# El nombre de la muestra se obtiene quitando el sufijo _Norm
sample_ids <- sub("_Norm$", "", norm_columns)

if (anyDuplicated(sample_ids)) {
  stop(
    "Hay identificadores de muestra duplicados después de ",
    "eliminar el sufijo '_Norm'. Revisa el Excel."
  )
}

message("Número de muestras _Norm: ", length(norm_columns))


# 6. COMPROBAR CORRESPONDENCIA ENTRE MUESTRAS Y METADATOS -------

metadata_ids <- metadata$patient_id

samples_without_metadata <- setdiff(sample_ids, metadata_ids)
metadata_without_samples <- setdiff(metadata_ids, sample_ids)

if (length(samples_without_metadata) > 0) {
  warning(
    "Muestras sin metadatos correspondientes: ",
    paste(samples_without_metadata, collapse = ", ")
  )
}

if (length(metadata_without_samples) > 0) {
  warning(
    "Pacientes de los metadatos sin columna _Norm: ",
    paste(metadata_without_samples, collapse = ", ")
  )
}

# Solo se analizan pacientes con muestra _Norm y metadatos
common_ids <- sample_ids[sample_ids %in% metadata_ids]

if (length(common_ids) < 3) {
  stop(
    "Hay menos de tres pacientes con coincidencia entre ",
    "las columnas _Norm y los metadatos."
  )
}

norm_columns <- paste0(common_ids, "_Norm")

# Conservar el orden de las muestras en la matriz
metadata_analysis <- metadata[
  match(common_ids, metadata$patient_id),
  metadata_columns,
  drop = FALSE
]

rownames(metadata_analysis) <- metadata_analysis$patient_id

# Verificar grupos clínicos
valid_status <- c("HEALTHY", "AFFECTED")

unexpected_status <- setdiff(
  unique(na.omit(metadata_analysis$clinical_status)),
  valid_status
)

if (length(unexpected_status) > 0) {
  stop(
    "Valores no esperados en clinical_status: ",
    paste(unexpected_status, collapse = ", "),
    ". Se esperan HEALTHY y AFFECTED."
  )
}

if (anyNA(metadata_analysis$clinical_status)) {
  stop("Hay clinical_status vacíos en los pacientes seleccionados.")
}

if (!all(valid_status %in% metadata_analysis$clinical_status)) {
  stop("Se necesitan muestras de ambos grupos: HEALTHY y AFFECTED.")
}

metadata_analysis$clinical_status <- factor(
  metadata_analysis$clinical_status,
  levels = c("HEALTHY", "AFFECTED")
)

write.csv(
  metadata_analysis,
  file.path(output_tables, "patient_metadata_used.csv"),
  row.names = FALSE,
  na = ""
)


# 7. PREPARAR IDENTIFICADORES DE REGION -------------------------

# Columnas de anotación que pueden existir en el Excel
annotation_candidates <- c(
  "#CHROM",
  "START",
  "END",
  "REGION ID",
  "GEN"
)

available_annotations <- intersect(
  annotation_candidates,
  names(raw_data)
)

if ("REGION ID" %in% names(raw_data)) {
  region_id <- as.character(raw_data[["REGION ID"]])
} else {
  region_id <- paste0("region_", seq_len(nrow(raw_data)))
}

empty_region_id <- is.na(region_id) | trimws(region_id) == ""
region_id[empty_region_id] <- paste0(
  "region_",
  which(empty_region_id)
)

# Evitar identificadores duplicados de región
region_id <- make.unique(region_id)

# Matriz numérica de regiones por pacientes
norm_df <- raw_data[, norm_columns, drop = FALSE]

norm_df[] <- lapply(norm_df, function(x) {
  suppressWarnings(as.numeric(as.character(x)))
})

norm_matrix <- as.matrix(norm_df)
storage.mode(norm_matrix) <- "numeric"

rownames(norm_matrix) <- region_id
colnames(norm_matrix) <- common_ids

# Exportar matriz original _Norm
norm_export <- data.frame(
  region_id = rownames(norm_matrix),
  norm_matrix,
  check.names = FALSE
)

write.csv(
  norm_export,
  file.path(output_tables, "norm_matrix.csv"),
  row.names = FALSE,
  na = ""
)

# Guardar las anotaciones de las regiones
region_annotations <- data.frame(
  region_id = region_id,
  raw_data[, available_annotations, drop = FALSE],
  check.names = FALSE
)

write.csv(
  region_annotations,
  file.path(output_tables, "region_annotations.csv"),
  row.names = FALSE,
  na = ""
)


# 8. CONTROL DE CALIDAD -----------------------------------------

missing_by_sample <- colSums(is.na(norm_matrix))
detected_by_sample <- colSums(!is.na(norm_matrix))
total_by_sample <- colSums(norm_matrix, na.rm = TRUE)

qc_table <- data.frame(
  patient_id = colnames(norm_matrix),
  clinical_status = metadata_analysis[
    colnames(norm_matrix),
    "clinical_status"
  ],
  n_regions = nrow(norm_matrix),
  n_detected_regions = detected_by_sample,
  n_missing_regions = missing_by_sample,
  fraction_missing = missing_by_sample / nrow(norm_matrix),
  total_norm = total_by_sample,
  median_norm = apply(norm_matrix, 2, median, na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(
  qc_table,
  file.path(output_tables, "sample_quality_control.csv"),
  row.names = FALSE,
  na = ""
)

# Resumen de filtrado
filter_summary <- data.frame(
  metric = c(
    "Total de regiones en Excel",
    "Total de pacientes con metadatos",
    "Total de regiones con al menos un valor observado",
    "Regiones con todos los valores observados",
    "Regiones sin ningún valor observado"
  ),
  value = c(
    nrow(norm_matrix),
    ncol(norm_matrix),
    sum(rowSums(!is.na(norm_matrix)) > 0),
    sum(rowSums(!is.na(norm_matrix)) == ncol(norm_matrix)),
    sum(rowSums(!is.na(norm_matrix)) == 0)
  )
)

write.csv(
  filter_summary,
  file.path(output_tables, "filter_summary.csv"),
  row.names = FALSE
)


# 9. GRAFICO: TOTAL _Norm POR PACIENTE --------------------------

p_totals <- ggplot(
  qc_table,
  aes(
    x = reorder(patient_id, total_norm),
    y = total_norm,
    fill = clinical_status
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Total de valores _Norm por paciente",
    x = "Patient ID",
    y = "Suma de valores _Norm",
    fill = "Clinical status"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(output_figures, "norm_totals.png"),
  p_totals,
  width = 9,
  height = 7,
  dpi = 300
)


# 10. GRAFICO: DISTRIBUCIONES POR PACIENTE ----------------------

distribution_data <- data.frame(
  patient_id = rep(colnames(norm_matrix), each = nrow(norm_matrix)),
  value = as.vector(norm_matrix),
  stringsAsFactors = FALSE
)

distribution_data$clinical_status <- metadata_analysis[
  distribution_data$patient_id,
  "clinical_status"
]

distribution_data <- distribution_data[
  is.finite(distribution_data$value),
  ,
  drop = FALSE
]

p_distribution <- ggplot(
  distribution_data,
  aes(x = value, colour = clinical_status)
) +
  geom_density(linewidth = 0.8, na.rm = TRUE) +
  facet_wrap(~ patient_id, scales = "free", ncol = 4) +
  labs(
    title = "Distribución de valores _Norm por paciente",
    x = "Valor _Norm",
    y = "Densidad",
    colour = "Clinical status"
  ) +
  theme_bw(base_size = 9)

ggsave(
  file.path(output_figures, "norm_distributions.png"),
  p_distribution,
  width = 13,
  height = max(7, ceiling(ncol(norm_matrix) / 4) * 2.2),
  dpi = 300
)


# 11. CORRELACION ENTRE PACIENTES -------------------------------

# Correlación de Pearson usando regiones compartidas observadas
cor_matrix <- cor(
  norm_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)

write.csv(
  data.frame(patient_id = rownames(cor_matrix), cor_matrix,
             check.names = FALSE),
  file.path(output_tables, "sample_correlation.csv"),
  row.names = FALSE,
  na = ""
)

annotation_col <- metadata_analysis[
  colnames(norm_matrix),
  c(
    "pathology",
    "gene_variant",
    "relationship",
    "genetic_status",
    "clinical_status"
  ),
  drop = FALSE
]

# Convertir anotaciones a factores para los colores del heatmap
annotation_col[] <- lapply(annotation_col, as.factor)

png(
  file.path(output_figures, "sample_correlation_heatmap.png"),
  width = 2400,
  height = 2200,
  res = 220
)

pheatmap::pheatmap(
  cor_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  main = "Correlación entre pacientes",
  border_color = NA,
  na_col = "grey85",
  fontsize = 8,
  angle_col = 45
)

dev.off()


# 12. PCA -------------------------------------------------------

# Para PCA se excluyen regiones sin valores o sin variación.
# Los valores ausentes restantes se imputan con la mediana de
# cada región exclusivamente para esta visualización exploratoria.

valid_rows_pca <- rowSums(!is.na(norm_matrix)) > 0

pca_matrix <- norm_matrix[valid_rows_pca, , drop = FALSE]

if (nrow(pca_matrix) < 2) {
  stop("No hay suficientes regiones para calcular PCA.")
}

row_medians <- apply(pca_matrix, 1, median, na.rm = TRUE)

for (i in seq_len(nrow(pca_matrix))) {
  missing_i <- is.na(pca_matrix[i, ])
  if (any(missing_i)) {
    pca_matrix[i, missing_i] <- row_medians[i]
  }
}

# Eliminar regiones sin variación
row_variance <- apply(pca_matrix, 1, var)
pca_matrix <- pca_matrix[
  is.finite(row_variance) & row_variance > 0,
  ,
  drop = FALSE
]

if (nrow(pca_matrix) < 2) {
  stop("No hay suficientes regiones variables para calcular PCA.")
}

# Transformación exploratoria
pca_input <- log2(pca_matrix + 1)

pca_result <- prcomp(
  t(pca_input),
  center = TRUE,
  scale. = TRUE
)

pca_coordinates <- data.frame(
  patient_id = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  stringsAsFactors = FALSE
)

pca_coordinates <- dplyr::left_join(
  pca_coordinates,
  metadata_analysis,
  by = "patient_id"
)

variance_explained <- 100 * (
  pca_result$sdev^2 / sum(pca_result$sdev^2)
)

write.csv(
  pca_coordinates,
  file.path(output_tables, "pca_coordinates.csv"),
  row.names = FALSE,
  na = ""
)

p_pca <- ggplot2::ggplot(
  pca_coordinates,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    colour = clinical_status,
    shape = genetic_status
  )
) +
  ggplot2::geom_point(size = 3) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = patient_id),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::labs(
    title = "PCA de valores _Norm",
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "%)"
    ),
    colour = "Clinical status",
    shape = "Genetic status"
  ) +
  ggplot2::theme_classic()

print(p_pca)

ggplot2::ggsave(
  filename = file.path(output_figures, "norm_pca1.png"),
  plot = p_pca,
  width = 10,
  height = 7,
  dpi = 300
)

pca_coordinates$pathology <- factor(
  pca_coordinates$pathology,
  levels = c("HEALTHY CONTROLS", "ALPS", "ALPS-LIKE")
)

pca_coordinates$clinical_status <- factor(
  pca_coordinates$clinical_status,
  levels = c("HEALTHY", "AFFECTED")
)

p_pca <- ggplot2::ggplot(
  pca_coordinates,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    colour = pathology,
    shape = clinical_status
  )
) +
  ggplot2::geom_point(size = 4) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = patient_id),
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "HEALTHY" = 16,
      "AFFECTED" = 17
    ),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "PCA de valores _Norm",
    subtitle = "Controles sanos, ALPS y ALPS-LIKE",
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "%)"
    ),
    colour = "Pathology",
    shape = "Clinical status"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position = "right"
  )

print(p_pca)

ggplot2::ggsave(
  filename = file.path(output_figures, "norm_pca2.png"),
  plot = p_pca,
  width = 12,
  height = 8,
  dpi = 300
)

# 13. ANALISIS DIFERENCIAL EXPLORATORIO CON LIMMA ---------------

# El modelo compara AFFECTED frente a HEALTHY.
# El análisis se limita a regiones con valores observados en
# todos los pacientes para no imputar datos en el modelo diferencial.
#
# IMPORTANTE:
# log2(_Norm + 1) es una transformación exploratoria. Si _Norm son
# conteos de citosinas metiladas sin denominador de cobertura,
# no equivale a analizar proporciones de metilación beta.

complete_rows <- rowSums(!is.na(norm_matrix)) == ncol(norm_matrix)

model_matrix <- norm_matrix[complete_rows, , drop = FALSE]

if (nrow(model_matrix) < 2) {
  stop(
    "Hay menos de dos regiones completas en todas las muestras. ",
    "Revisa los valores ausentes o define un filtro adecuado."
  )
}

if (any(model_matrix < 0, na.rm = TRUE)) {
  stop(
    "Se encontraron valores _Norm negativos. ",
    "Revisa los datos antes de usar log2(_Norm + 1)."
  )
}

model_log <- log2(model_matrix + 1)

group <- factor(
  metadata_analysis[colnames(model_log), "clinical_status"],
  levels = c("HEALTHY", "AFFECTED")
)

design <- model.matrix(~ group)
colnames(design) <- make.names(colnames(design))

if (!("groupAFFECTED" %in% colnames(design))) {
  stop("No se pudo construir el contraste AFFECTED frente a HEALTHY.")
}

fit <- limma::lmFit(model_log, design)
fit <- limma::eBayes(fit)

differential_results <- limma::topTable(
  fit,
  coef = "groupAFFECTED",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

differential_results$region_id <- rownames(differential_results)

# Añadir anotaciones de las regiones
region_annotations_unique <- region_annotations[
  match(differential_results$region_id, region_annotations$region_id),
  ,
  drop = FALSE
]

differential_results <- dplyr::left_join(
  differential_results,
  region_annotations_unique,
  by = "region_id"
)

# FDR de Benjamini-Hochberg calculado por limma
differential_results$significant_FDR <- (
  !is.na(differential_results$adj.P.Val) &
    differential_results$adj.P.Val < fdr_threshold
)

# Exportar todos los resultados
write.csv(
  differential_results,
  file.path(output_tables, "differential_norm_all_regions.csv"),
  row.names = FALSE,
  na = ""
)

# Exportar solo regiones significativas por FDR
significant_results <- differential_results |>
  dplyr::filter(significant_FDR) |>
  dplyr::arrange(adj.P.Val)

write.csv(
  significant_results,
  file.path(output_tables, "differential_norm_FDR05.csv"),
  row.names = FALSE,
  na = ""
)

message(
  "Regiones con FDR < ",
  fdr_threshold,
  ": ",
  nrow(significant_results)
)


# 14. VOLCANO PLOT: FDR -----------------------------------------

volcano_data <- differential_results |>
  dplyr::mutate(
    neg_log10_FDR = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
    significance = ifelse(
      significant_FDR,
      "FDR < 0.05",
      "Not significant"
    )
  )

# Etiquetar las regiones con FDR significativo; si no hay,
# mostrar las regiones con menor FDR como referencia.
labels_to_show <- volcano_data |>
  dplyr::filter(significant_FDR) |>
  dplyr::arrange(adj.P.Val) |>
  dplyr::slice_head(n = 15)

if (nrow(labels_to_show) == 0) {
  labels_to_show <- volcano_data |>
    dplyr::arrange(adj.P.Val) |>
    dplyr::slice_head(n = 10)
}

p_volcano <- ggplot(
  volcano_data,
  aes(
    x = logFC,
    y = neg_log10_FDR,
    colour = significance
  )
) +
  geom_point(alpha = 0.75, size = 1.7, na.rm = TRUE) +
  geom_hline(
    yintercept = -log10(fdr_threshold),
    linetype = "dashed"
  ) +
  ggrepel::geom_text_repel(
    data = labels_to_show,
    aes(label = region_id),
    size = 2.8,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  labs(
    title = "Análisis diferencial de valores _Norm",
    subtitle = "AFFECTED frente a HEALTHY; significación por FDR",
    x = "log2 fold change (limma)",
    y = "-log10(FDR)",
    colour = "Significance"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(output_figures, "volcano_norm.png"),
  p_volcano,
  width = 10,
  height = 7,
  dpi = 300
)


# 15. HEATMAP DE REGIONES VARIABLES / DIFERENCIALES -------------

# Si hay regiones significativas, usar las de menor FDR.
# Si ninguna alcanza FDR < 0.05, mostrar regiones variables
# como visualización exploratoria, sin afirmar significación.

if (nrow(significant_results) > 0) {
  
  heatmap_regions <- head(
    significant_results$region_id,
    top_regions_heatmap
  )
  
  heatmap_title <- paste0(
    "Regiones con FDR < ",
    fdr_threshold
  )
  
} else {
  
  row_variances <- apply(model_log, 1, var)
  row_variances <- row_variances[is.finite(row_variances)]
  
  heatmap_regions <- names(
    sort(row_variances, decreasing = TRUE)
  )[seq_len(min(top_regions_heatmap, length(row_variances)))]
  
  heatmap_title <- "Regiones más variables (sin regiones FDR significativas)"
}

heatmap_matrix <- model_log[
  intersect(heatmap_regions, rownames(model_log)),
  ,
  drop = FALSE
]

if (nrow(heatmap_matrix) >= 2) {
  
  # Estandarizar cada región para mostrar patrones relativos
  heatmap_z <- t(scale(t(heatmap_matrix)))
  heatmap_z[!is.finite(heatmap_z)] <- 0
  
  png(
    file.path(output_figures, "top_norm_regions_heatmap.png"),
    width = 2600,
    height = 2200,
    res = 220
  )
  
  pheatmap::pheatmap(
    heatmap_z,
    annotation_col = annotation_col,
    show_colnames = TRUE,
    show_rownames = nrow(heatmap_z) <= 60,
    labels_col = colnames(heatmap_z),
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    border_color = NA,
    main = heatmap_title,
    fontsize = 8,
    angle_col = 45
  )
  
  dev.off()
  
} else {
  warning("No hay suficientes regiones para generar el heatmap.")
}


# 16. RESUMEN FINAL ---------------------------------------------

summary_table <- data.frame(
  metric = c(
    "Samples analyzed",
    "HEALTHY samples",
    "AFFECTED samples",
    "Total regions in matrix",
    "Complete regions used for limma",
    "Regions with FDR < 0.05",
    "FDR threshold",
    "Differential model",
    "Input columns"
  ),
  value = c(
    ncol(norm_matrix),
    sum(metadata_analysis$clinical_status == "HEALTHY"),
    sum(metadata_analysis$clinical_status == "AFFECTED"),
    nrow(norm_matrix),
    nrow(model_matrix),
    nrow(significant_results),
    fdr_threshold,
    "limma on log2(_Norm + 1), AFFECTED vs HEALTHY",
    "Only columns ending in _Norm"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(output_tables, "analysis_summary.csv"),
  row.names = FALSE
)

message("==============================================")
message("ANALISIS FINALIZADO")
message("Tablas:  ", normalizePath(output_tables))
message("Figuras: ", normalizePath(output_figures))
message(
  "Regiones significativas por FDR < ",
  fdr_threshold,
  ": ",
  nrow(significant_results)
)
message("==============================================")