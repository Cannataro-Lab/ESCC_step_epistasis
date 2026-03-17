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

plot_data <- mutation_rates %>%
  mutate(is_highlight = gene %in% genes_of_interest) %>%
  # Filter to genes with non-zero rates to avoid log errors
  filter(normal_mu > 0, cancer_mu > 0)

text_size <- 15
geom_text_size <- text_size * (5/14)

mut_rate_plot_revised <- ggplot(plot_data, aes(x = normal_mu, y = cancer_mu)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.2, color = "grey70", size = 1) +
  geom_point(data = filter(plot_data, is_highlight), aes(color = gene), size = 4) +
  geom_text_repel(data = filter(plot_data, is_highlight), aes(label = gene, color = gene),
                  # fontface = "bold", 
                  box.padding = 0.5,
                  size = geom_text_size,
                  direction = "both",
                  force = 5,
                  force_pull = 1,
                  max.overlaps = Inf) +
  # Set limits slightly above zero
  scale_x_continuous(limits = c(0, 7.5e-7), labels = scientific) +
  scale_y_continuous(limits = c(0, 3.5e-6), labels = scientific) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Expected mutation burden in CHNE",
       y = "Expected mutation burden in ESCC") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = text_size), 
        axis.title = element_text(size = text_size), 
        panel.grid = element_line(color = scales::alpha("grey92", .4)),
        aspect.ratio = 1)


ggsave("output_data/fig_1_mutrates.png", mut_rate_plot_revised, width = 8, height = 6)
ggsave("output_data/fig_1_mutrates.pdf", mut_rate_plot_revised, width = 8, height = 6)
ggsave("output_data/fig_1_mutrates.jpg", mut_rate_plot_revised, width = 8, height = 6)




