required_packages <- c("dplyr", "readr", "Seurat")
invisible(lapply(required_packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) stop("Missing package: ", x)
}))

source("config.R")
dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(config$gse202379_seurat_rds)
required_meta <- c("cell.annotation", "Disease.status")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta)) stop("Missing metadata: ", paste(missing_meta, collapse = ", "))
if (!"SRC" %in% rownames(obj)) stop("SRC is absent from the RNA assay.")

obj$Group <- factor(
  ifelse(obj$Disease.status == "Healthy control", "Healthy", "MASLD"),
  levels = c("Healthy", "MASLD")
)
Idents(obj) <- "cell.annotation"

cell_counts <- as.data.frame(table(obj$cell.annotation, obj$Group))
names(cell_counts) <- c("cell_type", "group", "n_nuclei")
cell_props <- cell_counts |>
  dplyr::group_by(group) |>
  dplyr::mutate(proportion = n_nuclei / sum(n_nuclei)) |>
  dplyr::ungroup()
readr::write_csv(cell_props, file.path(config$results_dir, "GSE202379_cell_composition.csv"))

# Nucleus-level summaries are descriptive. Biological inference must use donor/sample IDs.
rna_data <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
src_by_cell_type <- lapply(unique(as.character(obj$cell.annotation)), function(cell_type) {
  cells_h <- colnames(obj)[obj$cell.annotation == cell_type & obj$Group == "Healthy"]
  cells_m <- colnames(obj)[obj$cell.annotation == cell_type & obj$Group == "MASLD"]
  if (!length(cells_h) || !length(cells_m)) return(NULL)
  mean_h <- mean(expm1(rna_data["SRC", cells_h]))
  mean_m <- mean(expm1(rna_data["SRC", cells_m]))
  data.frame(
    cell_type = cell_type,
    n_healthy = length(cells_h), n_masld = length(cells_m),
    mean_healthy = mean_h, mean_masld = mean_m,
    log2_fold_change = log2((mean_m + 1) / (mean_h + 1))
  )
})
src_by_cell_type <- dplyr::bind_rows(src_by_cell_type)
readr::write_csv(src_by_cell_type, file.path(config$results_dir, "GSE202379_SRC_descriptive_by_cell_type.csv"))
capture.output(sessionInfo(), file = file.path(config$results_dir, "GSE202379_expression_sessionInfo.txt"))

