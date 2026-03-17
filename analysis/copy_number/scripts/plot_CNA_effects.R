setwd("analysis/copy_number/")

library(data.table)
library(ggplot2)
library(scales)
library(patchwork)

# Define CNA events of interest based on previous literature
cna_of_interest <- c("FGF4", "CDKN2B", "SOX2", "MYC", "FOXA1", "GBE1",
                     "CCND1", "TEAD4", "EIF4EBP1", "EGFR", "FGFR3", "FZD6",
                     "PIK3CA", "FAM135B", "BCHE", "SDHA")
type_of_event <- c("increase", "decrease", "increase", "increase", "increase", "decrease",
                   "increase", "increase", "increase", "increase", "increase", "increase",
                   "increase", "increase", "increase", "increase")

interest_dt <- data.table(gene1 = cna_of_interest, event_type = type_of_event)

# Assemble results for plotting
# may need to adjust some things here if you set up your directory differently
dir_wt  <- file.path(ROOT, "untracked/outputs/")
dir_mut <- file.path(ROOT, "untracked/outputs/")

# List rds for all chromosomes
list_chr_rds <- function(base_dir) {
  list.files(base_dir,
             pattern = "ESCC_chr.*_iterative_rounds\\.rds$",
             recursive = TRUE,
             full.names = TRUE)
}

# Extract chr from filename (expects something like "..._chr21_...")
extract_chr <- function(path) {
  m <- regmatches(path, regexpr("_chr[0-9XY]+_", path))
  gsub("^_chr|_$", "", m)
}

# Extract round-1 single-term LRT p-value for a given gene/event
# Uses: round_results[[1]] which is the output of run_round() for round 1.
# We compress results to one row per (range_id,event_type) by taking the max loglik across runs
get_round1_stats <- function(res_obj) {
  if (is.null(res_obj$round_results) || length(res_obj$round_results) < 1) return(NULL)
  r1 <- as.data.table(res_obj$round_results[[1]])
  if (!nrow(r1)) return(NULL)
  
  # neutral base loglik stored in res_obj$base_loglik (if present)
  base_ll <- res_obj$base_loglik
  if (is.null(base_ll) || is.na(base_ll)) base_ll <- NA_real_
  
  # summarize to one row per candidate
  out <- r1[, .(
    gene1      = gene1[1],
    event_type = event_type[1],
    range_id   = range_id[1],
    chr        = chr[1],
    start      = start[1],
    end        = end[1],
    si         = median(si, na.rm = TRUE),
    loglik     = max(loglik, na.rm = TRUE),
    num_samples= max(num_samples, na.rm = TRUE),
    base_loglik= base_ll
  ), by = .(range_id, event_type)]
  
  # LRT vs base (df = 1) if base_ll exists
  out[, delta_ll := loglik - base_loglik]
  out[, lrt := 2 * pmax(delta_ll, 0)]
  out[, p_lrt := ifelse(is.finite(lrt) & !is.na(base_loglik),
                        pchisq(lrt, df = 1, lower.tail = FALSE),
                        NA_real_)]
  out[, neglog10p := -log10(pmax(p_lrt, .Machine$double.xmin))]
  out
}

# Determine if a gene/event was selected in the final iterative set
get_fixed_flags <- function(res_obj) {
  fx <- as.data.table(res_obj$fixed_events)
  if (!nrow(fx)) return(data.table(range_id=character(), event_type=character()))
  unique(fx[, .(range_id, event_type)])
}

# Get first round selected for a given gene/event by scanning fixed_events order
get_first_selected_round <- function(res_obj) {
  fx <- as.data.table(res_obj$fixed_events)
  if (!nrow(fx)) return(data.table(range_id=character(), event_type=character(), first_round=integer()))
  fx[, first_round := seq_len(.N)]
  fx[, .(first_round = min(first_round)), by = .(range_id, event_type)]
}

summarize_cohort <- function(cohort_name, base_dir) {
  files <- list_chr_rds(base_dir)
  if (!length(files)) stop("No *_iterative_rounds.rds found under: ", base_dir)
  
  all_interest_hits <- list()
  
  for (f in files) {
    chr <- extract_chr(f)
    res <- readRDS(f)
    
    r1 <- get_round1_stats(res)
    if (is.null(r1)) next
    
    fx <- get_fixed_flags(res)
    fr <- get_first_selected_round(res)
    
    # Join interest list by gene + expected event_type
    r1_interest <- r1[interest_dt, on = .(gene1, event_type), nomatch = 0]
    
    if (!nrow(r1_interest)) next
    
    r1_interest[, cohort := cohort_name]
    r1_interest[, chr_file := chr]
    
    # Mark if selected
    r1_interest[fx, is_selected := TRUE, on = c("range_id","event_type")]
    r1_interest[is.na(is_selected), is_selected := FALSE]
    
    # Add first round if selected
    r1_interest <- merge(r1_interest, fr, by = c("range_id","event_type"), all.x = TRUE)
    
    all_interest_hits[[length(all_interest_hits)+1]] <- r1_interest
  }
  
  out <- rbindlist(all_interest_hits, fill = TRUE)
  if (!nrow(out)) {
    warning("No matches for interest genes in cohort: ", cohort_name)
  }
  out
}

# Extract TP53 results
wt_tbl  <- summarize_cohort("TP53_WT",  dir_wt)
mut_tbl <- summarize_cohort("TP53_MUT", dir_mut)

both <- rbindlist(list(wt_tbl, mut_tbl), fill = TRUE)

# Ensure each (gene,event,cohort) is unique; if duplicates exist 
# keep the interval with strongest evidence (max neglog10p)
both_best <- both[order(-neglog10p)][
  , .SD[1],
  by = .(gene1, event_type, cohort)
]
both_best[, si := signif(si, 4)]
both_best[, delta_ll := signif(delta_ll, 4)]
both_best[, neglog10p := signif(neglog10p, 4)]

fwrite(both_best, "untracked/outputs/cna_interest_summary_table_tp53.tsv", sep = "\t")

# Extract NFE2L2 results
wt_tbl  <- summarize_cohort("NFE2L2_WT",  dir_wt)
mut_tbl <- summarize_cohort("NFE2L2_MUT", dir_mut)

both <- rbindlist(list(wt_tbl, mut_tbl), fill = TRUE)

# Ensure each (gene,event,cohort) is unique; if duplicates exist 
# keep the interval with strongest evidence (max neglog10p)
both_best <- both[order(-neglog10p)][
  , .SD[1],
  by = .(gene1, event_type, cohort)
]
both_best[, si := signif(si, 4)]
both_best[, delta_ll := signif(delta_ll, 4)]
both_best[, neglog10p := signif(neglog10p, 4)]

fwrite(both_best, "untracked/outputs/cna_interest_summary_table_nfe2l2.tsv", sep = "\t")


# Load CNA results if needed
tp53_cna_si_diff <- read_tsv("untracked/outputs/cna_interest_summary_table_tp53.tsv")
nfe2l2_cna_si_diff <- read_tsv("untracked/outputs/cna_interest_summary_table_nfe2l2.tsv")

tp53_cna_si_diff_wide <- tp53_cna_si_diff %>% 
  select(gene1, event_type, cohort, si, loglik, p_lrt) %>% 
  pivot_wider(names_from  = cohort,
              values_from = c(si, loglik, p_lrt),
              names_glue  = "{.value}_{cohort}")

tp53_cna_plot <- ggplot(tp53_cna_si_diff_wide) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey") +
  geom_point(aes(x = si_TP53_WT, y = si_TP53_MUT, shape = event_type, fill = event_type), size = 4) +
  geom_text_repel(aes(x = si_TP53_WT, y = si_TP53_MUT, label=gene1), box.padding = 0.5, min.segment.length = 2) +
  scale_shape_manual(values = c("increase" = 24, "decrease" = 25), guide = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = c("increase" = "darkblue", "decrease" = "darkorange"), guide = guide_legend(reverse = TRUE)) +
  expand_limits(x = 0, y = 0) +
  labs(x = "Selection for CNA in TP53 wildtype samples", y = "Selection for CNA in TP53 mutant samples",
       shape = "CNA event type",  fill = "CNA event type") +
  theme_bw()


nfe2l2_cna_si_diff_wide <- nfe2l2_cna_si_diff %>% 
  select(gene1, event_type, cohort, si, loglik, p_lrt) %>% 
  pivot_wider(names_from  = cohort,
              values_from = c(si, loglik, p_lrt),
              names_glue  = "{.value}_{cohort}")

nfe2l2_cna_plot <- ggplot(nfe2l2_cna_si_diff_wide) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey") +
  geom_point(aes(x = si_NFE2L2_WT, y = si_NFE2L2_MUT, shape = event_type, fill = event_type), size = 4) +
  geom_text_repel(aes(x = si_NFE2L2_WT, y = si_NFE2L2_MUT, label=gene1), box.padding = 0.3, min.segment.length = 2) +
  scale_shape_manual(values = c("increase" = 24, "decrease" = 25), guide = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = c("increase" = "darkblue", "decrease" = "darkorange"), guide = guide_legend(reverse = TRUE)) +  expand_limits(x = 0, y = 0) +
  labs(x = "Selection for CNA in NFE2L2 wildtype samples", y = "Selection for CNA in NFE2L2 mutant samples",
       shape = "CNA event type",  fill = "CNA event type") +
  theme_bw()



# We also jointly refitted the model for chromosomes where more than one driver 
# gene is on the same chromosome arm
joint_si_tbl <- read_csv("untracked/outputs/joint_chr3_chr8_chr11_selection_coefficients.tsv")
joint_override <- joint_si_tbl %>%
  transmute(
    gene1 = gene1.y,
    event_type = requested_event_type,
    cohort,
    chromosome,
    si_joint = si
  ) %>%
  left_join(
    loo_tbl %>%
      transmute(
        gene1 = dropped_gene,
        event_type = dropped_event_type,
        cohort,
        chromosome,
        delta_loglik_joint = delta_loglik,
        p_lrt_joint = p_lrt
      ),
    by = c("gene1", "event_type", "cohort", "chromosome")
  )

tp53_override <- tp53_cna_si_diff %>%
  left_join(
    joint_override,
    by = c("gene1", "event_type", "cohort")
  ) %>%
  mutate(
    si_for_plot = ifelse(!is.na(si_joint), si_joint, si),
    p_for_plot = ifelse(!is.na(p_lrt_joint), p_lrt_joint, p_lrt),
    delta_ll_for_plot = ifelse(!is.na(delta_loglik_joint), delta_loglik_joint, delta_ll),
    neglog10p_for_plot = -log10(pmax(p_for_plot, .Machine$double.xmin)),
    used_joint_model = !is.na(si_joint)
  )

nfe2l2_override <- nfe2l2_cna_si_diff %>%
  left_join(
    joint_override,
    by = c("gene1", "event_type", "cohort")
  ) %>%
  mutate(
    si_for_plot = ifelse(!is.na(si_joint), si_joint, si),
    p_for_plot = ifelse(!is.na(p_lrt_joint), p_lrt_joint, p_lrt),
    delta_ll_for_plot = ifelse(!is.na(delta_loglik_joint), delta_loglik_joint, delta_ll),
    neglog10p_for_plot = -log10(pmax(p_for_plot, .Machine$double.xmin)),
    used_joint_model = !is.na(si_joint)
  )

tp53_cna_si_diff_wide <- tp53_override %>% 
  select(gene1, event_type, cohort, si_for_plot, p_for_plot) %>% 
  pivot_wider(names_from  = cohort,
              values_from = c(si_for_plot, p_for_plot),
              names_glue  = "{.value}_{cohort}")

nfe2l2_cna_si_diff_wide <- nfe2l2_override %>% 
  select(gene1, event_type, cohort, si_for_plot, p_for_plot) %>% 
  pivot_wider(names_from  = cohort,
              values_from = c(si_for_plot, p_for_plot),
              names_glue  = "{.value}_{cohort}")

# Plotting for each stratified analysis for the joint refitting
text_size <- 15
geom_text_size <- text_size * (5/14)

tp53_cna_plot <- ggplot(tp53_cna_si_diff_wide) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey") +
  geom_point(aes(x = si_for_plot_TP53_WT, y = si_for_plot_TP53_MUT, shape = event_type, fill = event_type), size = 4) +
  geom_text_repel(aes(x = si_for_plot_TP53_WT, y = si_for_plot_TP53_MUT, label=gene1), 
                  box.padding = 0.3, min.segment.length = 1, size = geom_text_size) +
  scale_shape_manual(values = c("increase" = 24, "decrease" = 25), guide = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = c("increase" = "darkblue", "decrease" = "darkorange"), guide = guide_legend(reverse = TRUE)) +
  expand_limits(x = 0, y = 0) +
  labs(x = "Selection for CNA in TP53 wildtype samples", y = "Selection for CNA in TP53 mutant samples",
       shape = "CNA event type",  fill = "CNA event type") +
  scale_x_continuous(breaks = c(0,1,2,3)) +
  scale_y_continuous(breaks = c(0,1,2,3)) +
  theme_bw() +
  theme(axis.text = element_text(size = text_size), 
        axis.title = element_text(size = text_size),
        panel.grid.minor = element_blank(),
        legend.position = "none")

nfe2l2_cna_plot <- ggplot(nfe2l2_cna_si_diff_wide) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey") +
  geom_point(aes(x = si_for_plot_NFE2L2_WT, y = si_for_plot_NFE2L2_MUT, shape = event_type, fill = event_type), size = 4) +
  geom_text_repel(aes(x = si_for_plot_NFE2L2_WT, y = si_for_plot_NFE2L2_MUT, label=gene1), 
                  box.padding = 0.4, min.segment.length = 0.5, size = geom_text_size) +
  scale_shape_manual(values = c("increase" = 24, "decrease" = 25), guide = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = c("increase" = "darkblue", "decrease" = "darkorange"), guide = guide_legend(reverse = TRUE)) +
  expand_limits(x = 0, y = 0) +
  labs(x = "Selection for CNA in NFE2L2 wildtype samples", y = "Selection for CNA in NFE2L2 mutant samples",
       shape = "CNA event type",  fill = "CNA event type") +
  scale_x_continuous(breaks = c(0,1,2,3)) +
  scale_y_continuous(breaks = c(0,1,2,3)) +
  theme_bw() +
  theme(axis.text = element_text(size = text_size), 
        axis.title = element_text(size = text_size),
        legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size),
        panel.grid.minor = element_blank())

combined_cna_plot <- tp53_cna_plot + 
  nfe2l2_cna_plot + 
  plot_layout(ncol = 2)

ggsave("combined_CNA_plot.pdf", plot = combined_cna_plot, units = "in", height = 6, width = 12)