library(tidyverse)
library(data.table)

ucla_data_url <- "https://cbioportal-datahub.s3.amazonaws.com/escc_ucla_2014.tar.gz"
if(!file.exists("input_data/raw_data/escc_ucla_2014.tar.gz")){
  download.file(ucla_data_url, destfile = "input_data/raw_data/escc_ucla_2014.tar.gz") 
}
untar("input_data/raw_data/escc_ucla_2014.tar.gz", exdir = "input_data/raw_data/")

ucla_data <- read_tsv("input_data/raw_data/escc_ucla_2014/data_mutations.txt") 
ucla_data <- ucla_data %>%
  select(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, Sequence_Source) %>%
  mutate(Pre_or_Pri = "Pri")

data.table::fwrite(ucla_data, "input_data/ucla.maf", sep = "\t")
