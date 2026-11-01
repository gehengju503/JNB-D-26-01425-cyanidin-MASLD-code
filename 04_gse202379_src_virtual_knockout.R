required_packages <- c("clusterProfiler", "org.Hs.eg.db", "readr", "scTenifoldKnk", "Seurat")
invisible(lapply(required_packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) stop("Missing package: ", x)
}))

source("config.R")
dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

obj <- readRDS(config$gse202379_seurat_rds)
Idents(obj) <- "cell.annotation"
hep <- subset(obj, idents = "Hepatocytes")
hep <- Seurat::FindVariableFeatures(hep, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
genes_use <- intersect(unique(c("SRC", Seurat::VariableFeatures(hep))), rownames(hep))
cells_use <- sample(colnames(hep), size = min(2000L, ncol(hep)))
counts <- Seurat::GetAssayData(hep, assay = "RNA", layer = "counts")[genes_use, cells_use, drop = FALSE]
if (!"SRC" %in% rownames(counts)) stop("SRC is absent after feature selection.")

ko <- scTenifoldKnk::scTenifoldKnk(
  countMatrix = counts, gKO = "SRC", qc = TRUE,
  nc_nNet = 10, nc_nCells = 500, nc_nComp = 3, td_K = 3
)
saveRDS(ko, file.path(config$results_dir, "SRC_scTenifoldKnk_result.rds"))

ko_table <- as.data.frame(ko$diffRegulation)
if (!all(c("gene", "p.adj") %in% names(ko_table))) {
  stop("Unexpected result columns: ", paste(names(ko_table), collapse = ", "))
}
readr::write_csv(ko_table, file.path(config$results_dir, "SRC_virtual_knockout_diffRegulation.csv"))

sig_symbols <- unique(ko_table$gene[!is.na(ko_table$p.adj) & ko_table$p.adj < 0.05])
id_map <- clusterProfiler::bitr(
  sig_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db
)
go_bp <- clusterProfiler::enrichGO(
  gene = unique(id_map$ENTREZID), OrgDb = org.Hs.eg.db::org.Hs.eg.db,
  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
  pvalueCutoff = 0.05, qvalueCutoff = 0.20, readable = TRUE
)
readr::write_csv(as.data.frame(go_bp), file.path(config$results_dir, "SRC_virtual_knockout_GO_BP.csv"))
capture.output(sessionInfo(), file = file.path(config$results_dir, "SRC_virtual_knockout_sessionInfo.txt"))

