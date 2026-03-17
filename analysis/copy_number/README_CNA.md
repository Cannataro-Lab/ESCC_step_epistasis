# ESCC CNA analysis

## Recommended repository layout

``` text
analysis/copy_number/
├── scripts/
│   ├── prep_escc_cna.R
│   ├── prep_escc_TP53_stratified.R
│   ├── prep_escc_NFE2L2_stratified.R
│   ├── iterative_escc_bootstrap.R
│   └── plot_CNA_effects.R
├── prepped_data/
│   ├── TCGA_ASCAT3_hg38.txt.gz
│   ├── TCGA_TP53_mut_status.csv
│   ├── TCGA_NFE2L2_mut_status.csv
│   └── TCGA_patient_info.txt
├── untracked/
│   ├── tcga_esca_clinical.tsv
│   └── outputs/
└── slurm_logs/ (if needed for HPC)
```

## Software requirements

This pipeline depends on development versions of several packages.

-   R \>= 4.3 recommended
-   `ces.refset.hg38` from `Townsend-Lab-Yale/ces.refset.hg38@dev`
-   `cancereffectsizeR` from `Townsend-Lab-Yale/cancereffectsizeR@dev`
-   `BISCUT` from `jeff-mandell/BISCUT-py3`
-   CRAN packages used in the scripts: `data.table`, `readr`, `dplyr`, `ggplot2`, `ggrepel`, `pbapply`, `RColorBrewer`, `tidyr`, `patchwork`, `remotes`

The pipeline uses internal `cancereffectsizeR` functions via `:::` and depends on dev branches:

-   `ces.refset.hg38`: dev branch, commit 4bca8091f1d2bda0ed36a904ade28e8f2ec7eace
-   `cancereffectsizeR`: dev branch, commit ea7e199ec04af74135ddf0166b7b53f2ee855f7e
-   The scripts still use internal `cancereffectsizeR` functions via `:::` because those are the functions currently used in the analysis
-   Updates will be made to the scripts once CNA functions are released with a new `cancereffectsizeR` version

## Input data assumptions

The scripts assume the following inputs already exist:

-   ASCAT3 allele-specific segment calls for TCGA samples in hg38 coordinates
-   TCGA ESCA clinical metadata
-   patient sex metadata for expected chrX copy handling

The prep script filters to ESCC samples by matching `diagnoses.primary_diagnosis` to `Squamous cell carcinoma`.

## Step 1: prep ESCC CNA inputs

``` bash
Rscript scripts/prep_escc_cna.R \
  --root analysis/copy_number \
```

Main output:

-   `untracked/outputs/prepped_cna.rds`

For stratified analyses run `prep_escc_TP53_stratified.R` and/or `prep_escc_NFE2L2_stratified.R` instead

## Step 2: run iterative chromosome-level analysis

Run a single chromosome manually:

``` bash
Rscript scripts/iterative_escc_bootstrap.R 8 \
  --root analysis/copy_number \
  --prepped-rds analysis/copy_number/untracked/outputs/prepped_cna.rds \
  --output-root analysis/copy_number/untracked/outputs \
```

Or submit all chromosomes as an array job.

Main outputs per chromosome:

-   `ESCC_chr<CHR>_iterative_rounds.rds`
-   round-by-round PNG plots