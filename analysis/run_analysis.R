library(tidyverse)
library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg19)
library(R.utils)


# Set working directory
# setwd()


# Download ESCC-META data to input_data/raw_data/

# Download liftover files from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/ and https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/
download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz", destfile = "input_data/hg38ToHg19.over.chain.gz")
gunzip("input_data/hg38ToHg19.over.chain.gz")

download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz", destfile = "input_data/hg18ToHg19.over.chain.gz")
gunzip("input_data/hg18ToHg19.over.chain.gz")


# Read in modified input data ----
ucla_data <- read_tsv("input_data/ucla.maf")
yokoyama_data <- read_tsv("input_data/yokoyama.maf")
yuan_data <- read_tsv("input_data/yuan.maf")
martincorena_data <- read_tsv("input_data/martincorena.maf")
liu_data <- read_tsv("input_data/liu.maf")

# Read in data from ESCC-META dataset ----
# (excluded PMID32929369 data because it includes lung samples) 
ds1 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/ICGC/mutations_hg19.csv")
ds2 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID22877736/mutations_hg18.csv") 
ds3 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID25151357/mutations_hg18.csv") 
ds4 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID25839328/mutations_hg19.csv")
ds5 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID26873401/mutations_hg19.csv")
ds6 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID27058444/mutations_hg19.csv")
ds7 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID28548104/mutations_hg19.csv")
ds8 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID28608921/mutations_hg19.csv")
ds9 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID30012096/mutations_hg19.csv")
ds10 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID30975989/mutations_hg19.csv")
ds11 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID31289612/mutations_hg19.csv")
ds12 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID32398863/mutations_hg19.csv")
ds14 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID32974170/mutations_hg19.csv")
ds15 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID34263978/mutations_hg19.csv")
ds16 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID34285259/mutations_hg38.csv") 
ds17 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/PMID34413060/mutations_hg19.csv")
ds18 <- read_csv("~/../data/tumor_data/ESCA/ESCC-META/mutational_list/tcga_28052061/mutations_hg19.csv")

# Read in covered regions ----
covered_regions_ucla <- "targeted_regions/covered_regions_from_UCLA_Lin_S14.bed"
covered_regions_yoko <- "targeted_regions/yokoyama_covered_regions.bed"
covered_regions_yuan <- "targeted_regions/SureSelect_All_Exon_V5_S04380110_Covered_hg19.bed"
covered_regions_mart <- "targeted_regions/martincorena_covered_regions.bed"

# Prepare data and filter out likely germline/false-positive calls ----
ucla_data_maf <- preload_maf(maf = ucla_data, refset = "ces.refset.hg19", keep_extra_columns = "Pre_or_Pri")
ucla_data_maf <- ucla_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Pre_or_Pri = "Pri")

yokoyama_data_maf <- preload_maf(maf = yokoyama_data, refset = "ces.refset.hg19", keep_extra_columns = "Pre_or_Pri")
yokoyama_data_maf <- yokoyama_data_maf %>%
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  filter(variant_type == "snv")

yuan_data_maf <- preload_maf(maf = yuan_data, refset = "ces.refset.hg19", keep_extra_columns = "Pre_or_Pri")
yuan_data_maf <- yuan_data_maf %>%
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  filter(variant_type == "snv")

martincorena_data_maf <- preload_maf(maf = martincorena_data, refset = "ces.refset.hg19", keep_extra_columns = "Pre_or_Pri")
martincorena_data_maf <- martincorena_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Pre_or_Pri = "Pre")

liu_data_maf <- preload_maf(maf = liu_data, refset = "ces.refset.hg19", keep_extra_columns = "Pre_or_Pri")
liu_data_maf <- liu_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Pre_or_Pri = "Pri")

ds1_maf <- preload_maf(maf = ds1, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds2_maf <- preload_maf(maf = ds2, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT", chain_file = "input_data/hg18ToHg19.over.chain")
ds3_maf <- preload_maf(maf = ds3, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT", chain_file = "input_data/hg18ToHg19.over.chain")
ds4_maf <- preload_maf(maf = ds4, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds5_maf <- preload_maf(maf = ds5, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds6_maf <- preload_maf(maf = ds6, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds7_maf <- preload_maf(maf = ds7, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds8_maf <- preload_maf(maf = ds8, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds9_maf <- preload_maf(maf = ds9, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds10_maf <- preload_maf(maf = ds10, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds11_maf <- preload_maf(maf = ds11, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds12_maf <- preload_maf(maf = ds12, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds14_maf <- preload_maf(maf = ds14, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds15_maf <- preload_maf(maf = ds15, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds16_maf <- preload_maf(maf = ds16, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT", chain_file = "input_data/hg38ToHg19.over.chain")
ds17_maf <- preload_maf(maf = ds17, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")
ds18_maf <- preload_maf(maf = ds18, refset = "ces.refset.hg19", sample_col = "sample_id", chr_col = "CHROM", start_col = "POS", ref_col = "REF", tumor_allele_col = "ALT")

# Assign data from ESCC-META to primary tumor group ----
ds1_maf <- ds1_maf %>% mutate(Pre_or_Pri = "Pri")
ds2_maf <- ds2_maf %>% mutate(Pre_or_Pri = "Pri")
ds3_maf <- ds3_maf %>% mutate(Pre_or_Pri = "Pri")
ds4_maf <- ds4_maf %>% mutate(Pre_or_Pri = "Pri")
ds5_maf <- ds5_maf %>% mutate(Pre_or_Pri = "Pri")
ds6_maf <- ds6_maf %>% mutate(Pre_or_Pri = "Pri")
ds7_maf <- ds7_maf %>% mutate(Pre_or_Pri = "Pri")
ds8_maf <- ds8_maf %>% mutate(Pre_or_Pri = "Pri")
ds9_maf <- ds9_maf %>% mutate(Pre_or_Pri = "Pri")
ds10_maf <- ds10_maf %>% mutate(Pre_or_Pri = "Pri")
ds11_maf <- ds11_maf %>% mutate(Pre_or_Pri = "Pri")
ds12_maf <- ds12_maf %>% mutate(Pre_or_Pri = "Pri")
ds14_maf <- ds14_maf %>% mutate(Pre_or_Pri = "Pri")
ds15_maf <- ds15_maf %>% mutate(Pre_or_Pri = "Pri")
ds16_maf <- ds16_maf %>% mutate(Pre_or_Pri = "Pri")
ds17_maf <- ds17_maf %>% mutate(Pre_or_Pri = "Pri")
ds18_maf <- ds18_maf %>% mutate(Pre_or_Pri = "Pri")

# Specify WES and WGS samples from Zhang et al. ----
ds4_wgs_maf <- ds4_maf %>% filter(Unique_Patient_Identifier %in% 
                                    c("PMID25839328;1N01-VS-1T01",
                                      "PMID25839328;1N02-VS-1T02",
                                      "PMID25839328;1N03-VS-1T03",
                                      "PMID25839328;1N04-VS-1T04",
                                      "PMID25839328;1N05-VS-1T05",
                                      "PMID25839328;3N01-VS-3T01",
                                      "PMID25839328;3N02-VS-3T02",
                                      "PMID25839328;3N03-VS-3T03",
                                      "PMID25839328;3N04-VS-3T04",
                                      "PMID25839328;3N05-VS-3T05",
                                      "PMID25839328;3N06-VS-3T06",
                                      "PMID25839328;3N07-VS-3T07",
                                      "PMID25839328;3N08-VS-3T08",
                                      "PMID25839328;3N09-VS-3T09"))
ds4_wes_maf <- ds4_maf %>% filter(!Unique_Patient_Identifier %in% 
                                    c("PMID25839328;1N01-VS-1T01",
                                      "PMID25839328;1N02-VS-1T02",
                                      "PMID25839328;1N03-VS-1T03",
                                      "PMID25839328;1N04-VS-1T04",
                                      "PMID25839328;1N05-VS-1T05",
                                      "PMID25839328;3N01-VS-3T01",
                                      "PMID25839328;3N02-VS-3T02",
                                      "PMID25839328;3N03-VS-3T03",
                                      "PMID25839328;3N04-VS-3T04",
                                      "PMID25839328;3N05-VS-3T05",
                                      "PMID25839328;3N06-VS-3T06",
                                      "PMID25839328;3N07-VS-3T07",
                                      "PMID25839328;3N08-VS-3T08",
                                      "PMID25839328;3N09-VS-3T09"))

# Specify WES and WGS samples from Qin et al. ----
ds6_wes_maf <- ds6_maf %>% filter(grepl("QC", Unique_Patient_Identifier))
ds6_wgs_maf <- ds6_maf %>% filter(grepl("ASM", Unique_Patient_Identifier))

# Check for sample duplicates and sample contamination (already reviewed) ----
# mafs <- list(ucla = ucla_data_maf, 
#              yuan = yuan_data_maf,
#              yokoyama = yokoyama_data_maf,
#              martincorena = martincorena_data_maf,
#              liu = liu_data_maf,
#              ds1 = ds1_maf,
#              ds2 = ds2_maf,
#              ds3 = ds3_maf,
#              ds4 = ds4_maf,
#              ds5 = ds5_maf,
#              ds6 = ds6_maf,
#              ds7 = ds7_maf,
#              ds8 = ds8_maf,
#              ds9 = ds9_maf,
#              ds10 = ds10_maf,
#              ds11 = ds11_maf,
#              ds12 = ds12_maf,
#              ds14 = ds14_maf,
#              ds15 = ds15_maf,
#              ds16 = ds16_maf,
#              ds17 = ds17_maf,
#              ds18 = ds18_maf) 
# possible_dups <- check_sample_overlap(maf_list = mafs) 


# Set up cancereffectsizeR analysis ----
cesa <- CESAnalysis(refset = "ces.refset.hg19")


# Load in all maf files 
cesa <- load_maf(cesa, maf = ucla_data_maf, coverage = "targeted", sample_data_cols = "Pre_or_Pri", maf_name = "ucla",
                 covered_regions_name = "ucla", covered_regions = covered_regions_ucla, covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = yuan_data_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "yuan",
                 covered_regions_name = "yuan", covered_regions = covered_regions_yuan, covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = yokoyama_data_maf, coverage = "targeted", sample_data_cols = "Pre_or_Pri", maf_name = "yokoyama",
                 covered_regions_name = "yoko", covered_regions = covered_regions_yoko, covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = martincorena_data_maf, coverage = "targeted", sample_data_cols = "Pre_or_Pri", maf_name = "martincorena",
                 covered_regions_name = "mart", covered_regions = covered_regions_mart, covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = liu_data_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "liu")
cesa <- load_maf(cesa, maf = ds1_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds1")
cesa <- load_maf(cesa, maf = ds2_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds2")
cesa <- load_maf(cesa, maf = ds3_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds3")
cesa <- load_maf(cesa, maf = ds4_wgs_maf, coverage = "genome", sample_data_cols = "Pre_or_Pri", maf_name = "ds4")
cesa <- load_maf(cesa, maf = ds4_wes_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds4")
cesa <- load_maf(cesa, maf = ds5_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds5")
cesa <- load_maf(cesa, maf = ds6_wgs_maf, coverage = "genome", sample_data_cols = "Pre_or_Pri", maf_name = "ds6")
cesa <- load_maf(cesa, maf = ds6_wes_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds6")
cesa <- load_maf(cesa, maf = ds7_maf, coverage = "genome", sample_data_cols = "Pre_or_Pri", maf_name = "ds7")
cesa <- load_maf(cesa, maf = ds8_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds8")
cesa <- load_maf(cesa, maf = ds9_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds9")
cesa <- load_maf(cesa, maf = ds10_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds10")
cesa <- load_maf(cesa, maf = ds11_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds11")
cesa <- load_maf(cesa, maf = ds12_maf, coverage = "genome", sample_data_cols = "Pre_or_Pri", maf_name = "ds12")
cesa <- load_maf(cesa, maf = ds14_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds14")
cesa <- load_maf(cesa, maf = ds15_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds15")
cesa <- load_maf(cesa, maf = ds16_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds16")
cesa <- load_maf(cesa, maf = ds17_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds17")
cesa <- load_maf(cesa, maf = ds18_maf, coverage = "exome", sample_data_cols = "Pre_or_Pri", maf_name = "ds18")

save_cesa(cesa = cesa,file = "analysis/eso_cesa_before_generates.rds")


# Write supplementary table 1 outlining all samples used in the analysis ----
samples_used <- cesa@samples

samples_used <- samples_used %>%
  mutate(Study = case_when(maf_source == "ucla" ~ "Lin et al. (UCLA, cBioPortal)",
                           maf_source == "yuan" ~ "Yuan et al.",
                           maf_source == "yokoyama" ~ "Yokoyama et al.",
                           maf_source == "martincorena" ~ "Martincorena et al.",
                           maf_source == "liu" ~ "Liu et al.",
                           startsWith(maf_source, "ds") ~ "Li et al.")) %>%
  mutate(Sample = Unique_Patient_Identifier) %>%
  mutate(Normal_or_Tumor = case_when(Pre_or_Pri == "Pre" ~ "Normal",
                                     Pre_or_Pri == "Pri" ~ "Tumor")) %>%
  mutate(Notes = case_when(maf_source == "ds1" ~ "originally from ICGC",
                           maf_source == "ds2" ~ "originally published by Agrawal et al.(2012)",
                           maf_source == "ds3" ~ "originally published by Gao et al.(2014)",
                           maf_source == "ds4" ~ "originally published by Zhang et al.(2015)",
                           maf_source == "ds5" ~ "originally published by Sawada et al.(2016)",
                           maf_source == "ds6" ~ "originally published by Qin et al.(2016)",
                           maf_source == "ds7" ~ "originally published by Chang et al.(2017)",
                           maf_source == "ds8" ~ "originally published by Dai et al.(2017)",
                           maf_source == "ds9" ~ "originally published by Guo et al.(2018)",
                           maf_source == "ds10" ~ "originally published by Yan et al.(2019)",
                           maf_source == "ds11" ~ "originally published by Urabe et al.(2019)",
                           maf_source == "ds12" ~ "originally published by Cui et al.(2020)",
                           maf_source == "ds14" ~ "originally published by Mangalaparthi et al.(2018)",
                           maf_source == "ds15" ~ "originally published by Takemoto et al.(2021)",
                           maf_source == "ds16" ~ "originally published by Erkizan et al.(2021)",
                           maf_source == "ds17" ~ "originally published by Hirata et al.(2021)",
                           maf_source == "ds18" ~ "originally from TCGA",
                           !startsWith(maf_source, "ds") ~ "N/A")) %>%
  mutate(Coverage = coverage) %>%
  select(Study, Sample, Normal_or_Tumor, Coverage, Notes)

write_tsv(samples_used, "output_data/Supplementary_table_1.tsv")


# Create compound variant table ----
# Get consensus coverage across whichever samples you want to include.
# Here, we use all WES/TGS, but you could choose to exclude some if they don't cover the genes of interest well.
all_cov = c(cesa$coverage_ranges$exome, cesa$coverage_ranges$targeted)

# Exclude "exome", since typically "exome+" is what's applicable
all_cov = all_cov[! names(all_cov) == 'exome'] 
all_cov = Reduce(GenomicRanges::intersect, all_cov)

# Define however many genes of interest you want
genes = c("TP53", 
          "NOTCH1", 
          "NOTCH2", 
          "ERBB4", 
          "NFE2L2", 
          "PIK3CA", 
          "CDKN2A.p16INK4a",
          "CDKN2A.p14arf",
          "ARID1A", 
          "FAT1", 
          "EGFR",
          "ERBB2",
          "FBXW7",
          "FGFR3",
          "RB1",
          "SMAD4",
          "SOX2")

# gr argument will make only universally covered variants get returned
variants <- select_variants(cesa, genes = genes, gr = all_cov)

# Further filter variants table based on COSMIC oncogene/TSG classification (exclude nonrecurrent except nonsense for TSGs).
top_TP53 <- variants[gene == "TP53" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_NOTCH1 <- variants[gene == "NOTCH1" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_NOTCH2 <- variants[gene == "NOTCH2" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_ERBB4 <- variants[gene == "ERBB4" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_NFE2L2 <- variants[gene == "NFE2L2" & (maf_prevalence > 1)]
top_PIK3CA <- variants[gene == "PIK3CA" & (maf_prevalence >1)]
top_CDKN2A.p16INK4a <- variants[gene == "CDKN2A.p16INK4A" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_CDKN2A.p14arf <- variants[gene == "CDKN2A.p14arf" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_FAT1 <- variants[gene == "FAT1" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_EGFR <- variants[gene == "EGFR" & (maf_prevalence >1)]
top_ERBB2 <- variants[gene == "ERBB2" & (maf_prevalence >1)]
top_FBXW7 <- variants[gene == "FBXW7" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_FGFR3 <- variants[gene == "FGFR3" & (maf_prevalence >1)]
top_RB1 <- variants[gene == "RB1" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_SMAD4 <- variants[gene == "SMAD4" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_SOX2 <- variants[gene == "SOX2" & (maf_prevalence >1)]

for_comp <- rbind(top_TP53, 
                  top_NOTCH1, 
                  top_NOTCH2, 
                  top_ERBB4, 
                  top_NFE2L2, 
                  top_PIK3CA, 
                  top_CDKN2A.p16INK4a, 
                  top_CDKN2A.p14arf, 
                  top_FAT1, 
                  top_EGFR,
                  top_ERBB2,
                  top_FBXW7,
                  top_FGFR3,
                  top_RB1,
                  top_SMAD4,
                  top_SOX2)

# Filter out genes with maf prevalence less than 25 (looking at the data, 25 seemed like a sensible threshold to exclude non-significant results)
comp_genes_high_prevalence <- for_comp %>% 
  group_by(gene) %>% 
  summarize(occurrence = sum(maf_prevalence)) %>% 
  filter(occurrence > 25) 
for_comp <- for_comp %>%
  filter(gene %in% comp_genes_high_prevalence$gene)

# Define compound variants to find cancer effect sizes at the gene level and not for individual variants
compound <- define_compound_variants(cesa = cesa, variant_table = for_comp, by = "gene", merge_distance = Inf)


# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates ----
cesa <- gene_mutation_rates(cesa, covariates = "ESCA", samples = cesa$samples[Pre_or_Pri=="Pre"], save_all_dndscv_output = T)
cesa <- gene_mutation_rates(cesa, covariates = "ESCA", samples = cesa$samples[Pre_or_Pri=="Pri"], save_all_dndscv_output = T)



# Use hg19 reference set
RefCDS <- ces.refset.hg19$RefCDS

dndscv_gene_names <- cesa$gene_rates$gene

# Find number of synonymous sites
nsyn_sites <- sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

# Find number of samples in normal and tumor tissue
samples_in_normal <- length(unique(cesa$dNdScv_results$rate_grp_1$annotmuts$sampleID ))
samples_in_cancer <- length(unique(cesa$dNdScv_results$rate_grp_2$annotmuts$sampleID ))

mut_rate_df <- tibble(gene = cesa$dNdScv_results$rate_grp_1$genemuts$gene_name,
                      exp_normal_mu = cesa$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
                      exp_cancer_mu = cesa$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv)

mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

# Mu = ((expected mu)/(number of synonymous sites))/(samples)
mut_rate_df <- mut_rate_df %>% 
  mutate(normal_mu = (exp_normal_mu / n_syn_sites) / samples_in_normal) %>%
  mutate(cancer_mu = (exp_cancer_mu / n_syn_sites) / samples_in_cancer) %>%
  mutate(cancer_greater = cancer_mu > normal_mu) 

# Check that mutation rate in tumor tissue (stage 0->2) is greater than mutation rate in normal tissue (stage 0->1)
mut_rate_df$cancer_greater %>% table()

mut_rate_df <- mut_rate_df %>%
  mutate(mut_rate_normal = normal_mu) %>% # mutation rate from "stage 0->1"
  mutate(mut_rate_cancer = cancer_mu - normal_mu) # mutation rate from "stage 1->2"

mutation_rates <- mut_rate_df

# Get proportions of mutation rates ----
mut_rates_for_p <- mut_rate_df %>%
  select(gene,cancer_mu,normal_mu) %>% 
  mutate(p_1 = normal_mu / cancer_mu) %>% 
  mutate(p_2 = 1 - p_1)

set_cancer_rates <- mut_rate_df %>%
  select(gene, rate = cancer_mu) %>%
  data.table::setDT()

cesa <- clear_gene_rates(cesa = cesa)

# Now all samples have highest rates
cesa <- set_gene_rates(cesa = cesa, rates = set_cancer_rates, missing_genes_take_nearest = T) 


# Infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis ----
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "Eso-SCC")
cesa <- trinuc_mutation_rates(cesa, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)


# New sequential selection method ----
index_by_state = list()
name_by_state = list()
ordering_col = 'Pre_or_Pri'
ordering = c('Pre', 'Pri')

if(is.null(names(ordering))) {
  if (length(unlist(ordering)) == length(ordering)) {
    names(ordering) = unlist(ordering)
  } else {
    names(ordering) = 1:length(ordering)
  }
}

for (i in 1:length(ordering)) {
  for (j in 1:length(ordering[[i]])) {
    index_by_state[[ordering[[i]][j]]] = i
    name_by_state[[ordering[[i]][j]]] = names(ordering)[i]
  }
}

samples = cancereffectsizeR:::select_samples(cesa, samples=cesa$samples)
sample_index_table = samples[, .(Unique_Patient_Identifier = Unique_Patient_Identifier,
                                 group_index = unlist(index_by_state[samples[[ordering_col]]]), 
                                 group_name = unlist(name_by_state[samples[[ordering_col]]]))]

source("analysis/new_sequential_lik.R")

for(comp_ind in 1:length(compound)){
  this_comp <- compound[comp_ind, ]
  this_gene <- unlist(unique(this_comp$snv_info$genes))[1]
  these_props <- mut_rates_for_p[mut_rates_for_p$gene %in% this_gene, c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  
  if(length(this_gene) != 1){
    this_gene <- unlist(str_split(this_gene[1], "\\."))
    this_gene <- this_gene[1]
  }
  
  cat("Running gene:", this_gene, "\n")
  
  if (length(these_props) != 2 || any(is.na(these_props))) {
    stop(paste("Missing or invalid mutation rates for", this_gene))
  }
  
  cesa <- ces_variant(cesa = cesa, 
                      variants = this_comp, 
                      model = sequential_lik_dev, 
                      lik_args = list(sample_index = sample_index_table, 
                                      sequential_mut_prop = these_props), 
                      optimizer_args = list(method = 'L-BFGS-B', 
                                            lower = 1e-3, 
                                            upper = 1e9), # define optimizer arguments to avoid failed convergence 
                      # (cancereffectsizeR suppresses warnings about failed convergence)
                      return_fit = TRUE, # return fit for loglikelihood and confidence intervals 
                      run_name = this_gene,
                      conf = 0.95
  )
}

selection_results_step <- rbind(cesa@selection_results$TP53,
                                cesa@selection_results$NOTCH1,
                                cesa@selection_results$NOTCH2,
                                cesa@selection_results$NFE2L2,
                                cesa@selection_results$PIK3CA,
                                cesa@selection_results$FAT1,
                                cesa@selection_results$FBXW7,
                                cesa@selection_results$RB1)

# Calculate non-step-specific selection intensity for LRT----

# use mutation rates consistent with the method above
mu_zero_to_one <- mutation_rates %>%
  select(gene, rate = normal_mu) %>%
  data.table::setDT()
mu_zero_to_two <- mutation_rates %>%
  select(gene, rate = cancer_mu) %>% # this is just the sum of \mu_{0\to1} and \mu_{1\to2}
  data.table::setDT()


cesa <- clear_gene_rates(cesa = cesa)

cesa <- set_gene_rates(cesa = cesa, rates = mu_zero_to_one, missing_genes_take_nearest = T, samples = cesa$samples[Pre_or_Pri=="Pre"]) 
cesa <- set_gene_rates(cesa = cesa, rates = mu_zero_to_two, missing_genes_take_nearest = T, samples = cesa$samples[Pre_or_Pri=="Pri"]) 

# calculate selection intensities using the simpler one-step model
cesa <- ces_variant(cesa = cesa, 
                    variants = compound, 
                    return_fit = TRUE,
                    run_name = "simple_model",
                    conf = 0.95)


# Calculate likelihood ratio ----
genes <- c("TP53", "NOTCH1", "NOTCH2", "NFE2L2", "PIK3CA", "FAT1", "FBXW7", "RB1")

loglik_step <- selection_results_step$loglikelihood

loglik_simple <- cesa@selection_results$simple_model$loglikelihood

loglik_df <- data.frame(
  gene = genes,
  loglik_step = loglik_step,
  loglik_simple = loglik_simple)

loglik_df <- loglik_df %>%
  mutate(
    loglik_step = ifelse(loglik_step > 1e5, NA, loglik_step), # just making sure that the loglikelihoods are realistic
    loglik_simple = ifelse(loglik_simple > 1e5, NA, loglik_simple), # not some large positive number (happens when convergence failed)
    LRT_stat = -2 * (loglik_simple - loglik_step),
    p_value = pchisq(LRT_stat, df = 1, lower.tail = FALSE), 
    p_less_0.5 = ifelse(p_value < 0.05, TRUE, FALSE))


# Clear gene rates and calculate gene rates for all samples (not separated by normal and tumor) for epistasis ----
cesa <- clear_gene_rates(cesa)
cesa <- gene_mutation_rates(cesa, covariates = "ESCA", save_all_dndscv_output = T)

dndscv_gene_names <- cesa$gene_rates$gene
nsyn_sites <- sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

samples_in_all <- length(unique(cesa$dNdScv_results$rate_grp_1$annotmuts$sampleID ))

mut_rate_df <- tibble(gene = cesa$dNdScv_results$rate_grp_1$genemuts$gene_name,
                      exp_mu = cesa$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv)

mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df <- mut_rate_df %>% 
  mutate(total_mu = (exp_mu / n_syn_sites) / samples_in_all) %>%
  select(gene, total_mu) %>%
  data.table::setDT()

cesa <- clear_gene_rates(cesa = cesa)

mut_rate_df <- mut_rate_df %>%
  select(gene, rate = total_mu)
cesa <- set_gene_rates(cesa = cesa, rates = mut_rate_df, missing_genes_take_nearest = T) 


cesa <- ces_epistasis(cesa, variants = compound, run_name = "epistasis_compound_variants_all_samples")


save_cesa(cesa = cesa, file = "analysis/eso_cesa_after_analysis.rds")

