# E-drive recovery audit for GSE202379 and scTenifold code

## Recovered source locations

The E-drive search found the analysis project under:

`E:/后续实验规划/网络药理学文章/C3G-C3R-网络药理学/网络药理学-虚拟敲除过表达/C3G-C3R虚拟敲除过表达/`

The useful code was not saved as ordinary `.R` scripts. It was recovered from:

- `.Rhistory`
- `.Rproj.user/D4021377/sources/session-1883ab9a/5EB1C66C-contents`
- `.Rproj.user/D4021377/sources/session-1f212c19/2467B455-contents`
- `.Rproj.user/D4021377/ctx/*/history`

The public GSE202379 Seurat object was found beside the project as one extracted RDS and two approximately 2.84 GB gzip copies. A separate `GSE202379-sirt1-2.R` script was also found at `E:/后续实验规划/`, but it is a SIRT1/CytoTRACE/pseudotime workflow rather than the manuscript's SRC workflow.

## Scientific-integrity problems in the recovered cache

The cache is not suitable for direct public deposition. It contains:

1. Randomly simulated `mock_z` and `mock_p` values when scTenifold failed.
2. A manually constructed `go_data_exact` table rather than enrichment calculated from the actual significant genes.
3. Hard-coded annotations and figure labels that can survive even when calculations fail.
4. A pathway-scoring branch that adds one count to a list of downstream target genes rather than only to SRC.
5. An artificial `+0.8` addition to the overexpressed target score.
6. Repeated trial-and-error blocks, local F/G-drive paths, runtime package installation, and silent error handling.
7. Nucleus-level t tests without donor/sample-level replication control.

## Cleaning decisions

Scripts `03`-`05` retain only the documented analysis idea and recoverable parameters:

- hepatocytes;
- top 1,000 variable genes with SRC forced into the feature set;
- random seed 42;
- KO: up to 2,000 nuclei, `nc_nNet=10`, `nc_nCells=500`, `nc_nComp=3`, `td_K=3`;
- OE: SRC count `+1`, `nc_nNet=10`, `nc_nCells=200`, `nc_nComp=3`, `td_K=3`;
- adjusted-P threshold below 0.05;
- BH-corrected GO biological-process enrichment;
- Hallmark module scores with no artificial offsets.

No simulated results or silent fallback datasets are generated. A calculation failure stops execution.

## Mandatory next step

Regenerate Figure 9E-F and Figure 10 from the cleaned scripts and compare them with the submitted panels. If the regenerated results differ materially, update the figures, Results, legends, Supplement, response letter, and interpretation before resubmission. Donor-level metadata should be recovered from GEO and used for inferential differential-expression analyses; nucleus-level summaries should remain descriptive.
