library(cancereffectsizeR)
library(ggplot2)
library(cowplot)
library(data.table)


signature_table_twostage <- cesa$mutational_signatures$biological_weights
signature_table_twostage <- signature_table_twostage %>%
  select(!c("total_snvs", "sig_extraction_snvs", "group_avg_blended"))

tumor_pre <- cesa$samples %>%
  filter(Pre_or_Pri == "Pre") %>%
  select(Unique_Patient_Identifier)
tumor_pri <- cesa$samples %>%
  filter(Pre_or_Pri == "Pri") %>%
  select(Unique_Patient_Identifier)

# Check which signatures don't have median weights of zero
signature_table_twostage[tumor_pre, on = "Unique_Patient_Identifier"] %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(), median)) %>% 
  pivot_longer(cols = 1:18) %>% 
  arrange(desc(value)) %>%
  filter(value != 0)

signature_table_twostage[tumor_pri, on = "Unique_Patient_Identifier"] %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(), median)) %>% 
  pivot_longer(cols = 1:18) %>% 
  arrange(desc(value))%>%
  filter(value != 0)


# Remove signatures with all-zero weight
is_all_zero <- lapply(signature_table_twostage, function(x) all(x == 0))
cols_to_remove <- names(is_all_zero)[is_all_zero == TRUE]
signature_table_twostage <- signature_table_twostage[, -c(cols_to_remove), with = F]


# all tumors will be present in table, so no NAs
# not much point of looking at the normal tissue signatures since they were assumed same as primary
tumor_pri_signature <- signature_table_twostage[tumor_pri, on = "Unique_Patient_Identifier"]
tumor_pri_signature <- tumor_pri_signature %>%
  select(Unique_Patient_Identifier, SBS1, SBS2, SBS5, SBS13, SBS18)

longer <- tumor_pri_signature %>%
  pivot_longer(!Unique_Patient_Identifier, names_to = "signatures", values_to = "signature_weight")
longer$signatures <- as.factor(longer$signatures) %>%
  fct_relevel(c("SBS1","SBS2","SBS5","SBS13","SBS18"))


tumor_pri_signature_plot <- ggplot(data=longer, aes(x=signatures, y=signature_weight, fill=signatures)) +
  geom_violin(scale="width") +
  theme(plot.title = element_text(hjust = 0), 
        axis.text.x = element_text (hjust = 0.5, angle = 0),
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 18),
        title = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("COSMIC signature") +
  ylab("Signature weight") + ylim(0, 1) + scale_x_discrete() +
  ggtitle("Tumor") +
  theme(axis.line = element_line(color="black", linewidth = 0.5)) 
tumor_pri_signature_plot


tumor_pre_signature <- signature_table_twostage[tumor_pre, on = "Unique_Patient_Identifier"]
tumor_pre_signature <- tumor_pre_signature %>%
  select(Unique_Patient_Identifier, SBS1, SBS2, SBS5, SBS13, SBS18)

longer <- tumor_pre_signature %>%
  pivot_longer(!Unique_Patient_Identifier, names_to = "signatures", values_to = "signature_weight")
longer$signatures <- as.factor(longer$signatures) %>%
  fct_relevel(c("SBS1","SBS2","SBS5","SBS13","SBS18")) # only signatures without 0 median weight
  #fct_relevel(c("SBS1","SBS2","SBS5","SBS9","SBS10c","SBS10d","SBS13","SBS15","SBS17a","SBS17b","SBS18","SBS22","SBS36","SBS86","SBS87","SBS90","SBS93"))

tumor_pre_signature_plot <- ggplot(data=longer, aes(x=signatures, y=signature_weight, fill=signatures)) +
  geom_violin(scale="width") +
  theme(plot.title = element_text(hjust = 0), 
        axis.text.x = element_text (hjust = 0.5, angle = 0),
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 18),
        title = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("COSMIC signature") +
  ylab("Signature weight") + ylim(0, 1) + scale_x_discrete() +
  ggtitle("Normal") +
  theme(axis.line = element_line(color="black", linewidth = 0.5))
tumor_pre_signature_plot

signature_plots <- plot_grid(tumor_pre_signature_plot,tumor_pri_signature_plot, labels = c("", ""), ncol = 2)
signature_plots_with_label <- plot_grid(signature_plots, labels = c("B"), label_size = 24)

ggsave("output_data/fig_1b_signature_plots.png", signature_plots_with_label, width=15, height=9)



