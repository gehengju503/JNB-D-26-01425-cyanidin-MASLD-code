required_packages <- c("clusterProfiler", "dplyr", "readr", "scTenifoldNet", "Seurat")
invisible(lapply(required_packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) stop("Missing package: ", x)
}))

source("config.R")
dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

obj <- readRDS(config$gse202379_seurat_rds)
Idents(obj) <- "cell.annotation"
hep <- subset(obj, idents = "Hepatocytes")
hep$Group <- factor(
  ifelse(hep$Disease.status == "Healthy control", "Healthy", "MASLD"),
  levels = c("Healthy", "MASLD")
)
hep <- Seurat::FindVariableFeatures(hep, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
genes_use <- intersect(unique(c("SRC", Seurat::VariableFeatures(hep))), rownames(hep))
counts_wt <- Seurat::GetAssayData(hep, assay = "RNA", layer = "counts")[genes_use, , drop = FALSE]
if (!"SRC" %in% rownames(counts_wt)) stop("SRC is absent after feature selection.")

# Computational scenario reported in the manuscript: add one count to SRC in every hepatocyte.
counts_oe <- counts_wt
counts_oe["SRC", ] <- counts_oe["SRC", ] + 1
oe <- scTenifoldNet::scTenifoldNet(
  X = counts_wt, Y = counts_oe, qc = TRUE,
  nc_nNet = 10, nc_nCells = 200, nc_nComp = 3, td_K = 3
)
saveRDS(oe, file.path(config$results_dir, "SRC_scTenifoldNet_overexpression_result.rds"))
oe_table <- as.data.frame(oe$diffRegulation)
readr::write_csv(oe_table, file.path(config$results_dir, "SRC_virtual_overexpression_diffRegulation.csv"))

# Hallmark module scores. No artificial score offsets or simulated fallback values are used.
gmt <- clusterProfiler::read.gmt(config$hallmark_gmt)
hallmark <- split(gmt$gene, gmt$term)
hep_wt <- Seurat::AddModuleScore(hep, features = hallmark, name = "Hallmark_WT", seed = 42)

hep_oe <- hep
full_counts <- Seurat::GetAssayData(hep_oe, assay = "RNA", layer = "counts")
full_counts["SRC", ] <- full_counts["SRC", ] + 1
hep_oe <- Seurat::SetAssayData(hep_oe, assay = "RNA", layer = "counts", new.data = full_counts)
hep_oe <- Seurat::NormalizeData(hep_oe, assay = "RNA", verbose = FALSE)
hep_oe <- Seurat::AddModuleScore(hep_oe, features = hallmark, name = "Hallmark_OE", seed = 42)

score_summary <- function(meta, prefix, state) {
  score_cols <- grep(paste0("^", prefix, "[0-9]+$"), names(meta), value = TRUE)
  dplyr::bind_rows(lapply(seq_along(score_cols), function(i) {
    h <- meta[[score_cols[i]]][meta$Group == "Healthy"]
    m <- meta[[score_cols[i]]][meta$Group == "MASLD"]
    data.frame(
      state = state, pathway = names(hallmark)[i],
      mean_healthy = mean(h, na.rm = TRUE), mean_masld = mean(m, na.rm = TRUE),
      t_statistic_descriptive = unname(stats::t.test(m, h)$statistic)
    )
  }))
}

scores <- dplyr::bind_rows(
  score_summary(hep_wt@meta.data, "Hallmark_WT", "WT"),
  score_summary(hep_oe@meta.data, "Hallmark_OE", "SRC_plus_1")
)
readr::write_csv(scores, file.path(config$results_dir, "SRC_virtual_overexpression_hallmark_scores.csv"))
capture.output(sessionInfo(), file = file.path(config$results_dir, "SRC_virtual_overexpression_sessionInfo.txt"))
