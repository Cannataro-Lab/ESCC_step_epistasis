library(readxl)
library(tidyr)
library(dplyr)
library(data.table)

# Download Yuan et al. dataset from https://academic.oup.com/carcin/article/40/12/1445/5579375?login=true (Supplementary table 3)

yuan_file <- "input_data/raw_data/yuan_carcinogenesis_2019/bgz162_suppl_supplementary_table_3.xlsx"

yuan_data_1 <- read_xlsx(yuan_file, sheet=1) %>%
  mutate(Patient_No = 1)  %>% 
  pivot_longer(cols = N1:T4, names_to = "Normal_or_Tumor", values_to = "present") %>%
  select(!starts_with("LN")) %>%
  filter(present == 1)
yuan_data_2 <- read_xlsx(yuan_file, sheet=2) %>%
  mutate(Patient_No = 2)  %>% 
  pivot_longer(cols = N1:T6, names_to = "Normal_or_Tumor", values_to = "present") %>%
  select(!starts_with("LN")) %>%
  filter(present == 1)
yuan_data_3 <- read_xlsx(yuan_file, sheet=3) %>%
  mutate(Patient_No = 3)  %>% 
  pivot_longer(cols = N1:T4, names_to = "Normal_or_Tumor", values_to = "present") %>%
  select(!starts_with("LN")) %>%
  filter(present == 1)
yuan_data_4 <- read_xlsx(yuan_file, sheet=4) %>%
  mutate(Patient_No = 4)  %>% 
  pivot_longer(cols = N1:T4, names_to = "Normal_or_Tumor", values_to = "present") %>%
  select(!starts_with("LN")) %>%
  filter(present == 1)
yuan_data_5 <- read_xlsx(yuan_file, sheet=5) %>%
  mutate(Patient_No = 5)  %>% 
  pivot_longer(cols = N1:T6, names_to = "Normal_or_Tumor", values_to = "present") %>%
  select(!starts_with("LN")) %>%
  filter(present == 1)
yuan_data_6 <- read_xlsx(yuan_file, sheet=6) %>%
  mutate(Patient_No = 6)  %>% 
  pivot_longer(cols = N1:T4, names_to = "Normal_or_Tumor", values_to = "present") %>%
  select(!starts_with("LN")) %>%
  filter(present == 1)
yuan_data_7 <- read_xlsx(yuan_file, sheet=7) %>%
  mutate(Patient_No = 7)  %>% 
  pivot_longer(cols = N1:T4, names_to = "Normal_or_Tumor", values_to = "present") %>%
  select(!starts_with("LN")) %>%
  filter(present == 1)
yuan_data_8 <- read_xlsx(yuan_file, sheet=8) %>%
  mutate(Patient_No = 8)  %>% 
  pivot_longer(cols = N1:T4, names_to = "Normal_or_Tumor", values_to = "present") %>%
  select(!starts_with("LN")) %>%
  filter(present == 1)
yuan_data_9 <- read_xlsx(yuan_file, sheet=9) %>%
  mutate(Patient_No = 9)  %>% 
  pivot_longer(cols = N1:T4, names_to = "Normal_or_Tumor", values_to = "present") %>%
  select(!starts_with("LN")) %>%
  filter(present == 1)
yuan_data_10 <- read_xlsx(yuan_file, sheet=10) %>%
  mutate(Patient_No = 10)  %>% 
  pivot_longer(cols = N1:T4, names_to = "Normal_or_Tumor", values_to = "present") %>%
  select(!starts_with("LN")) %>%
  filter(present == 1)

all_yuan_data <- rbind(yuan_data_1, yuan_data_2, yuan_data_3, yuan_data_4, yuan_data_5, yuan_data_6, yuan_data_7, yuan_data_8, yuan_data_9, yuan_data_10)

all_yuan_data <- all_yuan_data %>%
  filter(Function =="exonic;splicing" | Function =="missense SNV" | Function =="stopgain" | Function =="stoploss" | Function =="synonymous SNV") %>% #filter for SNVs
  mutate(ref = substr(alter.,1,1)) %>%
  mutate(mut = substr(alter., 3,3)) %>%
  filter(ref != "-" & mut != "-") %>%
  mutate(Pre_or_Pri = ifelse(startsWith(Normal_or_Tumor, "N"), "Pre", "Pri")) %>% #assign "Pre" or "Pri" to normal and tumor samples
  mutate(Tumor_Sample_Barcode = ifelse(Pre_or_Pri=="Pre", paste0("P",Patient_No,"_",Normal_or_Tumor), paste0("P",Patient_No,"Tumor"))) %>% 
  select(Tumor_Sample_Barcode, Chromosome = Chr, Start_Position = Start, Reference_Allele = ref, Tumor_Seq_Allele2 = mut, Pre_or_Pri)

yuan_data <- all_yuan_data %>% # samples with severe dysplasia or carcinoma in situ share clones identical with tumors, so they were filtered out of normal samples
  filter(Tumor_Sample_Barcode != "P3_N1" & #SD
           Tumor_Sample_Barcode != "P3_N2" & #CIS
           Tumor_Sample_Barcode != "P3_N4" & #CIS
           Tumor_Sample_Barcode != "P6_N2" & #acanthosis (TP53 mutation identical with the primary tumor)
           Tumor_Sample_Barcode != "P6_N4" & #CIS
           Tumor_Sample_Barcode != "P6_N5" & #CIS
           Tumor_Sample_Barcode != "P9_N4") #CIS

data.table::fwrite(yuan_data, file = "input_data/yuan.maf", sep = "\t")

