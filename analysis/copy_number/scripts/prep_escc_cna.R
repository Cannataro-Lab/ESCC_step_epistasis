setwd("analysis/copy_number/")

library(devtools)
library(data.table)
library(remotes)
library(tidyverse)
library(readr)

remotes::install_github("Townsend-Lab-Yale/ces.refset.hg38@dev")
library(ces.refset.hg38)

# The version of ces.refset.hg38 that you need resides on the dev branch.
# Install with remotes::install_github("Townsend-Lab-Yale/ces.refset.hg38@dev")
stopifnot(packageVersion('ces.refset.hg38') >= as.package_version('1.3.0.9000'))

# cancereffectsizeR v3 also from dev branch. Rather than installing with install_github,
# clone the repo and use devtools as shown above.
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@dev")
library(cancereffectsizeR)
stopifnot(packageVersion('cancereffectsizeR') >= as.package_version('3.0.0.9000'))

# You also need a custom installation of BISCUT.
if(! require('BISCUT') || packageVersion('BISCUT') < as.package_version('1.2')) {
  remotes::install_github('jeff-mandell/BISCUT-py3')
}
library(BISCUT)


# Set directories
data_dir = 'prepped_data/' 
output_dir = 'untracked/outputs/'
dir.create(output_dir #, showWarnings = FALSE
)

# Read all TCGA cna calls
cna_calls <- fread(file.path(data_dir, 'TCGA_ASCAT3_hg38.txt.gz'))

# TCGA esca patients
tcga_esca_clinical <- read_tsv("untracked/tcga_esca_clinical.tsv")
tcga_escc_clinical <- tcga_esca_clinical %>%
  mutate(patient = cases.submitter_id) %>%
  filter(diagnoses.primary_diagnosis %like% 'Squamous cell carcinoma')

escc_cna_calls <- cna_calls %>%
  filter(patient_id %in% tcga_escc_clinical$cases.submitter_id)

# Load in sex (needed for expected chrX copy)
patient_info <- fread("prepped_data/TCGA_patient_info.txt")
escc_cna_calls[patient_info, sex := sex, on = 'patient_id']

# Filter to just patients in use (>99%)
patient_info <- patient_info[patient_id %in% escc_cna_calls$patient_id]
escc_cna_calls <- escc_cna_calls %>%
  filter(!is.na(sex)) %>%
  mutate(sex = case_when(sex == 'female' ~ 'F',
                         sex == 'male' ~ 'M'))
prepped_calls <- cancereffectsizeR:::prep_ASCAT3_segments(segments = escc_cna_calls, refset = 'ces.refset.hg38')

# Adjust cores as necessary for your system.
library(reticulate)
# reticulate::install_python(version = "3.11.11")
# py_config()

prepped_calls <- cancereffectsizeR:::call_large_events(prepped_calls = prepped_calls, arm_chr_threshold = .99, 
                                                       cores = 4, account_biscut_regions = TRUE,
                                                       run_biscut = FALSE)

# Run MutationalPatterns (takes a bit)
# Note that CN signatures were derived from hg38 data
cn_sig_def_file <- system.file('extdata/COSMIC_v3.4_CN_GRCh37.txt', package ='ces.refset.hg38')

signature_output <- cancereffectsizeR:::cn_signature_extraction(sig_def = cn_sig_def_file, 
                                                                cna_segments = prepped_calls$calls)

# For every sample and segment size class, get expected burden (under CN signature reconstruction)
# and rel_burden (proportion of overall burden among all samples within the size class).
# Output is 2-length list: separate data.table of burdens for increase_burden and decrease_burden
cna_burdens <- cancereffectsizeR:::cna_class_relative_rates(cna_calls = prepped_calls$calls,
                                                            cna_recon = signature_output$reconstruction,
                                                            ploidy_calls = prepped_calls$ploidy_calls)

# Save outputs for downstream analyses
# The unneeded plots that get included in MP output are cumulatively huge, so we'll remove
signature_output$mp_out$sim_decay_fig = 'deleted'
to_save = list(prepped_calls = prepped_calls, signature_output = signature_output,
               cna_burdens = cna_burdens)
saveRDS(to_save, file.path(output_dir, 'prepped_cna.rds'))