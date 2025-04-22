library(tidyverse)
library(data.table)
library(openxlsx)

# Download Martincorena et al. dataset from https://pmc.ncbi.nlm.nih.gov/articles/instance/6298579/bin/NIHMS80426-supplement-Supplementary_Table_2.xlsx
# Save file to input_data/raw_data/

martincorena_unmerged_table <- read.xlsx("input_data/raw_data/NIHMS80426-supplement-Supplementary_Table_2.xlsx", 
                                         sheet = 1, 
                                         startRow = 17)
martincorena_merged_table <- read.xlsx("input_data/raw_data/NIHMS80426-supplement-Supplementary_Table_2.xlsx", 
                                       sheet = 2, 
                                       startRow = 18)

martincorena_merged_tidy_samples <- martincorena_merged_table %>% 
  select(donor, gene, ntchange, merged) %>% # select columns relevant to analysis
  mutate(merged_index = 1:nrow(.)) %>%
  separate_rows(merged, sep = ";") %>% # break samples into multiple rows (tidy)
  mutate(merged = stringr::word(string = merged, start = 1, sep = ",")) # just keep the sample name (first "word")

martincorena_merged_tidy_samples <- martincorena_merged_tidy_samples %>%
  filter(!is.na(merged)) %>% 
  filter(merged != "too_far") %>%
  mutate(clone = paste(gene,ntchange))

merged_two_plus <- martincorena_merged_tidy_samples %>%
  count(merged) %>%
  filter(n > 1 ) %>%
  pull(merged)

martincorena_data <- martincorena_unmerged_table %>% 
  filter(!sampleID %in% merged_two_plus)

martincorena_data <- martincorena_data %>%
  select(Tumor_Sample_Barcode = sampleID, Chromosome = chr, Start_Position = pos, Reference_Allele = ref, Tumor_Seq_Allele2 = mut) %>%
  mutate(Pre_or_Pri = "Pre") #clean up data frame for CES analysis

data.table::fwrite(martincorena_data, 'input_data/martincorena.maf', sep = "\t")
