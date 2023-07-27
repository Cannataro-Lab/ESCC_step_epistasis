library(tidyverse)
library(data.table)
library(openxlsx)

yokoyama_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0811-x/MediaObjects/41586_2018_811_MOESM3_ESM.xlsx"
yokoyama_data <- read.xlsx(yokoyama_url, sheet = 8, startRow = 2)

yokoyama_data <- rename(yokoyama_data, c(type_of_mutations = Type.of.mutations,
                                         ref_allele = Reference.allele, 
                                         mut_allele = Mutant.allele))

yokoyama_data <- yokoyama_data %>%
  mutate(Chr = substr(Chr, 4, 5)) %>%
  filter(type_of_mutations == "nonsynonymous SNV" | type_of_mutations == "stopgain SNV" | type_of_mutations == "stoploss SNV" 
         | type_of_mutations == "synonymous SNV" | type_of_mutations == "splicing") #filter for SNVs

pre_yokoyama_data <- yokoyama_data %>%
  filter(Histology == "PNE") %>% #filter for physiologically normal samples
  select(Tumor_Sample_Barcode = Sample, Chromosome = Chr, Start_Position = Start, Reference_Allele = ref_allele, Tumor_Seq_Allele2 = mut_allele) %>%
  mutate(Pre_or_Pri = "Pre") 
pri_yokoyama_data <- yokoyama_data %>%
  filter(Histology == "ESCC") %>% #filter for ESCC samples
  select(Tumor_Sample_Barcode = Sample, Chromosome = Chr, Start_Position = Start, Reference_Allele = ref_allele, Tumor_Seq_Allele2 = mut_allele) %>%
  mutate(Pre_or_Pri = "Pri") 

all_yokoyama_data <- rbind(pre_yokoyama_data, pri_yokoyama_data) #combine both yokoyama data frames

#due to possible duplicates in the dataset, we are removing some problematic samples
samples_to_drop <- c("UPN37_1_13", "UPN37_1_14", "UPN37_1_15", "UPN37_1_17", "UPN37_1_8", "UPN37_1_12", "UPN37_1_16", "UPN37_1_18", "UPN37_2_3", 
                     "UPN37_2_10", "UPN37_2_1", "UPN49_1_6", "UPN58_1_19", "UPN58_1_16", "UPN58_1_17", "UPN58_1_14", "UPN36_1_8", "UPN36_1_2", 
                     "UPN36_1_19", "UPN36_1_3", "UPN55_1_2", "UPN55_1_6", "UPN55_1_14", "UPN55_1_5", "UPN55_1_24", "UPN55_1_21", "UPN55_1_3", 
                     "UPN55_1_22", "UPN55_1_7", "UPN55_1_17", "UPN55_1_18", "UPN55_1_9")

all_yokoyama_data <- all_yokoyama_data %>% 
  filter(!Tumor_Sample_Barcode %in% samples_to_drop)

data.table::fwrite(all_yokoyama_data, "input_data/yokoyama.maf", sep = "\t")
