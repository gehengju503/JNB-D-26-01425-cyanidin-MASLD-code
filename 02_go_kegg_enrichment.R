required_packages <- c("clusterProfiler", "enrichplot", "ggplot2", "org.Hs.eg.db", "readr")
invisible(lapply(required_packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) stop("Missing package: ", x)
}))

source("config.R")
dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)

gene_input <- readr::read_csv(config$gene_list_file, show_col_types = FALSE)
if (!"SYMBOL" %in% names(gene_input)) stop("Gene-list file must contain a SYMBOL column.")
symbols <- unique(stats::na.omit(trimws(gene_input$SYMBOL)))

id_map <- clusterProfiler::bitr(
  symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db
)
readr::write_csv(id_map, file.path(config$results_dir, "gene_id_mapping.csv"))

run_go <- function(ontology) {
  clusterProfiler::enrichGO(
    gene = unique(id_map$ENTREZID),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
}

go_results <- list(BP = run_go("BP"), CC = run_go("CC"), MF = run_go("MF"))
for (ontology in names(go_results)) {
  result_table <- as.data.frame(go_results[[ontology]])
  readr::write_csv(result_table, file.path(config$results_dir, paste0("GO_", ontology, ".csv")))
  if (nrow(result_table)) {
    p <- enrichplot::dotplot(go_results[[ontology]], showCategory = 15) +
      ggplot2::ggtitle(paste("GO", ontology))
    ggplot2::ggsave(
      file.path(config$results_dir, paste0("GO_", ontology, "_dotplot.png")),
      p, width = 9, height = 6, dpi = 600
    )
  }
}

kegg <- clusterProfiler::enrichKEGG(
  gene = unique(id_map$ENTREZID),
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  minGSSize = 10,
  maxGSSize = 500,
  use_internal_data = FALSE
)
readr::write_csv(as.data.frame(kegg), file.path(config$results_dir, "KEGG.csv"))
capture.output(sessionInfo(), file = file.path(config$results_dir, "enrichment_sessionInfo.txt"))

