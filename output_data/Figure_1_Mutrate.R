library(cancereffectsizeR)
library(data.table)
library(scales)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(tidyverse)

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}

# Load in cesa if necessary
# cesa <- load_cesa("analysis/eso_cesa_after_analysis.rds") 

genes_of_interest <- c("TP53", 
                       "NOTCH1", 
                       "NOTCH2", 
                       "NFE2L2", 
                       "PIK3CA", 
                       "FAT1", 
                       "FBXW7",
                       "RB1") #define genes of interest
mut_rates_specific_genes <- mutation_rates %>%
  mutate(highlight = ifelse(gene %in% genes_of_interest, TRUE, FALSE)) #highlight genes of interest
selected_mut_rates <- mutation_rates %>% 
  filter(gene %in% genes_of_interest)

selected_mut_rates_longer <- selected_mut_rates %>%
  pivot_longer(cols = c("normal_mu", "cancer_mu"), names_to = "progression", values_to = "mutation_rate")
selected_mut_rates_longer$progression <- factor(selected_mut_rates_longer$progression, levels = unique(selected_mut_rates_longer$progression))
only_cancer_rates <- selected_mut_rates_longer %>%
  filter(progression == "cancer_mu")

change_mut_rates_plot <- ggplot() + 
  geom_point(data = selected_mut_rates_longer, mapping = aes(x = progression, y = mutation_rate, group = gene, color = gene)) +
  geom_line(data = selected_mut_rates_longer, mapping = aes(x = progression, y = mutation_rate, group = gene, color = gene)) +
  geom_text_repel(data = only_cancer_rates, mapping = aes(x = 2.05, y = mutation_rate, label = gene, color = gene), hjust = -1, direction = "y", size = 5) +
  scale_y_continuous(labels=scientific) + 
  scale_x_discrete(labels=c("Conception to \nadult normal tissue", "Conception \nto tumor tissue")) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Trajectory", y = "Mutation rate", color = "Gene") +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        axis.title.x = element_text(vjust = -1.5), 
        legend.position = "none")

ggsave("output_data/fig_1_mutrates.png", change_mut_rates_plot, width = 8, height = 6)




