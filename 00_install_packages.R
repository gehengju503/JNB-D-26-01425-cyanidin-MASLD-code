cran_packages <- c(
  "broom", "dplyr", "ggplot2", "openxlsx", "readr", "readxl",
  "rms", "survey", "svyVGAM", "tibble", "ggrepel", "Seurat"
)

missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_packages <- c("clusterProfiler", "enrichplot", "org.Hs.eg.db")
missing_bioc <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE)

# scTenifoldKnk and scTenifoldNet are required for scripts 04 and 05.
# Install validated versions separately and record the exact source/commit.
