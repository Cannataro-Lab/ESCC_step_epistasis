library(tidyverse)
library(data.table)
library(openxlsx)

# Download Liu et al. dataset from https://www.gastrojournal.org/cms/10.1053/j.gastro.2017.03.033/attachment/37b9d627-f726-495b-8143-71c5f0e73f40/mmc2.xlsx
# Save file to input_data/raw_data/

liu_data <- read.xlsx("input_data/raw_data/mmc2.xlsx", sheet = 12, startRow = 2)

liu_data <- liu_data %>%
  mutate(type = sub(".*_", "", Sample)) %>%
  filter(type == "T") %>%
  select(Tumor_Sample_Barcode = Sample, Chromosome, Start_Position = Position, Reference_Allele = "Nucleotide.change.(ref.alt)", Tumor_Seq_Allele2 = X6) %>%
  mutate(Pre_or_Pri = "Pri") #clean up data frame for CES analysis

data.table::fwrite(liu_data, "input_data/liu.maf", sep = "\t")
