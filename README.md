# JNB-D-26-01425 reproducibility code

This repository accompanies the revised manuscript:

> Dietary cyanidin intake is inversely associated with metabolic dysfunction-associated steatotic liver disease, and cyanidin glycosides attenuate hepatocyte lipotoxicity

## What is included

- `R/01_nhanes_analysis.R`: survey-weighted descriptive analysis, primary and BMI-adjusted logistic models, restricted cubic splines, and analysis sample counts.
- `R/02_go_kegg_enrichment.R`: GO/KEGG over-representation analysis with Benjamini-Hochberg correction.
- `R/03_gse202379_src_expression.R`: GSE202379 cell-composition and descriptive cell-type SRC summaries.
- `R/04_gse202379_src_virtual_knockout.R`: recovered and cleaned SRC virtual-knockout workflow.
- `R/05_gse202379_src_virtual_overexpression.R`: recovered and cleaned SRC additive-overexpression and Hallmark scoring workflow.
- `config.example.R`: portable input/output paths.
- `CODE_RESPONSE_MATRIX.csv`: mapping from reviewer requests to available code and unresolved items.
- `UPLOAD_GUIDE_CN.md`: detailed Chinese instructions for GitHub and Zenodo deposition and manuscript updates.

## What is not included

Individual-level NHANES data, intervention RNA-seq analysis code, molecular-docking files, Cytoscape sessions, and MR scripts were not found as reusable source code in the reviewed folders. GSE202379/scTenifold source fragments were recovered from RStudio history/cache and rewritten without simulated fallback results. All recovered modules still require clean-environment rerunning and result reconciliation. See `CODE_RESPONSE_MATRIX.csv` and `RECOVERY_AUDIT.md`.

## Data inputs

1. NHANES 2017-2018 public-use files: <https://wwwn.cdc.gov/nchs/nhanes/>
2. FNDDS flavonoid data used to derive two-day cyanidin intake.
3. A locally prepared analysis table with the variable names documented in `data/README.md`.
4. A one-column gene-symbol file for enrichment analysis.

The local reviewed analysis table contained 3,882 participants; non-missing sample sizes were CAP/LSM 3,875, HSI 3,647, and FLI 1,756. These counts should be reconciled with the manuscript and participant flow diagram.

## Run

1. Install R 4.3 or later.
2. Copy `config.example.R` to `config.R` and edit the paths.
3. Install the packages listed in `R/00_install_packages.R` once.
4. Run `R/01_nhanes_analysis.R`.
5. Run `R/02_go_kegg_enrichment.R` if the gene-list input is available.
6. Run scripts `03`-`05` only after installing validated scTenifold packages and obtaining the public GSE202379 Seurat object and Hallmark GMT file.

Scripts do not install packages at run time and do not contain machine-specific absolute paths.

## Important validation note

The source folders contained analysis fragments in Word files named “错误代码-需校正”. The scripts here correct obvious programming defects and make the workflow portable, but the numerical outputs must be regenerated and compared with every manuscript table and figure before public release. In particular, confirm the final MASLD operational definition, exclusions, survey weights, factor coding, spline knots, and covariate definitions against the approved statistical analysis plan.

The recovered RStudio cache contained simulated fallback volcano data, a hard-coded GO table, modification of downstream target genes instead of SRC in one pathway-scoring branch, and an artificial `+0.8` score offset. None of those operations is included in the cleaned scripts. Previously generated Figure 9/10 panels must therefore be regenerated before submission.

## Suggested citation

After creating a Zenodo release, replace this section with the generated DOI and archive citation.
