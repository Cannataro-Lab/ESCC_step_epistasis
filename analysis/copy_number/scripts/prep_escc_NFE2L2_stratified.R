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

# Load the NFE2L2 status data
nfe2l2_status <- read_csv("prepped_data/TCGA_NFE2L2_mut_status.csv")

nfe2l2_mut_ids <- nfe2l2_status %>% filter(snv_status == "MUT") %>% pull(sample)
nfe2l2_wt_ids  <- nfe2l2_status %>% filter(snv_status == "WT") %>% pull(sample)

# Filter the calls into two separate objects
calls_mut <- escc_cna_calls %>% filter(patient_id %in% nfe2l2_mut_ids)
calls_wt  <- escc_cna_calls %>% filter(patient_id %in% nfe2l2_wt_ids)

# Define a function to process the cohorts
process_cna_cohort <- function(segment_data, label) {
  message("Processing cohort: ", label)
  prepped <- cancereffectsizeR:::prep_ASCAT3_segments(segments = segment_data, 
                                                      refset = 'ces.refset.hg38')
  
  prepped <- cancereffectsizeR:::call_large_events(prepped_calls = prepped, 
                                                   arm_chr_threshold = .99, 
                                                   cores = 4, 
                                                   account_biscut_regions = TRUE,
                                                   run_biscut = FALSE)
  
  cn_sig_def_file <- system.file('extdata/COSMIC_v3.4_CN_GRCh37.txt', 
                                 package ='ces.refset.hg38')
  
  sig_out <- cancereffectsizeR:::cn_signature_extraction(sig_def = cn_sig_def_file, 
                                                         cna_segments = prepped$calls)
  
  burdens <- cancereffectsizeR:::cna_class_relative_rates(cna_calls = prepped$calls,
                                                          cna_recon = sig_out$reconstruction,
                                                          ploidy_calls = prepped$ploidy_calls)
  
  sig_out$mp_out$sim_decay_fig = 'deleted'
  return(list(prepped_calls = prepped, signature_output = sig_out, cna_burdens = burdens))
}

# Run for NFE2L2 Mutant
res_mut <- process_cna_cohort(calls_mut, "NFE2L2 Mutant")
saveRDS(res_mut, file.path(output_dir, 'prepped_cna_NFE2L2_mut.rds'))

# Run for NFE2L2 Wild-type
res_wt <- process_cna_cohort(calls_wt, "NFE2L2 WT")
saveRDS(res_wt, file.path(output_dir, 'prepped_cna_NFE2L2_wt.rds'))
