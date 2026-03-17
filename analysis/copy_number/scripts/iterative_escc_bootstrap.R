#!/usr/bin/env Rscript

## ESCC iterative CNA analysis, one chromosome per run
## with BOOTSTRAP override for choosing other/noncancer vs driver

## Usage:
##   Rscript iterative_escc_bootstrap.R <CHR>
## Example:
##   Rscript iterative_escc_bootstrap.R 8


suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
  library(ggrepel)
  library(scales)
  library(pbapply)
  library(cancereffectsizeR)
  library(ces.refset.hg38)
})

ROOT <- "analysis/copy_number/"
setwd(ROOT)

## Stable temp dir (avoid node /tmp cleanup)
Sys.setenv(TMPDIR = file.path(ROOT, "tmp_R"))
dir.create(Sys.getenv("TMPDIR"), recursive = TRUE, showWarnings = FALSE)

## Parse chromosome argument 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript iterative_escc_bootstrap.R <CHR>\n",
       "Example: Rscript iterative_escc_bootstrap.R 8\n",
       "Allowed: 1-22 or X")
}
chr_id <- args[1]
proj   <- "ESCC"

message("Running ESCC all-gene iterative CNA analysis on chromosome ", chr_id)

## Optional: deal with scna_warning in cancereffectsizeR dev version
ns_env <- asNamespace("cancereffectsizeR")
unlockBinding <- get("unlockBinding", envir = asNamespace("base"))
unlockBinding("scna_warning", ns_env)
assign("scna_warning", value = function() {
  msg <- paste0(
    "Functions related to analysis of SCNAs (somatic copy number alterations) ",
    "are in beta. The underlying methodologies are still in development, ",
    "and documentation is incomplete. Please proceed with caution."
  )
  warning(msg, call. = FALSE)
  invisible(NULL)
}, envir = ns_env)
lockBinding("scna_warning", ns_env)

## Input: ESCC prepped calls 
proj_rds <- file.path(ROOT, "untracked/outputs/prepped_cna.rds")
if (!file.exists(proj_rds)) stop("Prepped CNA file not found at: ", proj_rds)

message("Loading prepped data from: ", proj_rds)
prepped_data <- readRDS(proj_rds)

## Robustly extract expected objects from the RDS
if (!is.null(prepped_data$prepped_calls) &&
    all(c("calls", "effective_coordinates") %in% names(prepped_data$prepped_calls))) {
  prepped_calls <- prepped_data$prepped_calls
} else if (all(c("calls", "effective_coordinates") %in% names(prepped_data))) {
  prepped_calls <- prepped_data
} else {
  stop("Unexpected structure in prepped_cna.rds.\n",
       "Expected either $prepped_calls$calls+$prepped_calls$effective_coordinates or top-level calls+effective_coordinates.\n",
       "Top-level names are: ", paste(names(prepped_data), collapse = ", "))
}

if (is.null(prepped_data$cna_burdens)) stop("prepped_cna.rds is missing `cna_burdens`.")

cna_calls        <- data.table::copy(prepped_calls$calls)
cna_burdens      <- prepped_data$cna_burdens
chr_coordinates  <- data.table::copy(prepped_calls$effective_coordinates)
chr_coordinates[, chr_label := paste0("chr", chr)]

## Output roots
output_root     <- file.path(ROOT, "untracked/outputs")
proj_output_dir <- file.path(output_root, proj)
dir.create(proj_output_dir, recursive = TRUE, showWarnings = FALSE)

## Segmentation rates for ESCC 
seg_rates <- cancereffectsizeR:::get_segmentation_rates(
  cna_calls   = cna_calls,
  cna_burdens = cna_burdens
)

## Gene-based disjoint intervals (ALL genes from refset) 
gc <- cancereffectsizeR:::get_disjoint_gene_coord(
  gene_coord = ces.refset.hg38$cancer_gene_coord
)
gc <- gc[, .(chr, start, end, gene1 = gene, gene1_anno = cancer_anno)]

chr_intervals <- gc[chr == chr_id]
if (!nrow(chr_intervals)) stop("No intervals found for chr ", chr_id, " in gene intervals.")
message("Number of gene intervals on chr", chr_id, ": ", nrow(chr_intervals))

## CNA rates for ESCC + chromosome
chr_rates <- cancereffectsizeR:::get_cna_rates(
  chr            = chr_id,
  cna_calls      = cna_calls,
  seg_rates      = seg_rates,
  rate_intervals = chr_intervals
)
chr_rates <- cancereffectsizeR:::clean_rates(chr_rates)
if (!nrow(chr_rates$rates)) stop("No rates left for ESCC chr ", chr_id, " after clean_rates().")
message("Number of rows in chr_rates$rates after clean_rates(): ", nrow(chr_rates$rates))

# Helper to run iterative rounds
run_round <- function(always_selected = data.table(),
                      rates          = NULL,
                      cores          = 1,
                      debug          = FALSE) {
  if (!is.data.table(always_selected)) stop("always_selected should be data.table.")
  for_events <- unique(rates$rates[, .(event_type, range_id)])
  
  if (always_selected[, .N] > 0) {
    stopifnot(all(c("range_id", "event_type") %in% names(always_selected)))
    for_events <- for_events[!always_selected, on = c("range_id", "event_type")]
    events_to_test <- split(always_selected$range_id, always_selected$event_type)
  } else {
    events_to_test <- list()
  }
  
  events_to_cycle <- split(for_events$range_id, for_events$event_type)
  
  output <- lapply(
    names(events_to_cycle),
    function(curr_event_type) {
      message("Testing models with one more ", curr_event_type, " effect term...")
      pbapply::pblapply(
        events_to_cycle[[curr_event_type]],
        function(curr_event) {
          curr_to_test <- copy(events_to_test)
          curr_to_test[[curr_event_type]] <- c(curr_to_test[[curr_event_type]], curr_event)
          cancereffectsizeR:::run_seg_multi(
            events_to_test = curr_to_test,
            rates          = rates,
            debug          = debug
          )
        },
        cl = cores
      )
    }
  )
  
  output <- unlist(output, recursive = FALSE, use.names = FALSE)
  
  si_out      <- rbindlist(lapply(output, "[[", "si"), idcol = "run")
  ll          <- sapply(output, function(x) x$info$loglik)
  num_samples <- sapply(output, function(x) length(x$info$included_samples))
  
  si_out[, loglik      := ll[run]]
  si_out[, num_samples := num_samples[run]]
  
  if (always_selected[, .N] > 0) {
    si_out[always_selected, is_fixed := TRUE, on = c("range_id", "event_type")]
    si_out[is.na(is_fixed), is_fixed := FALSE]
  } else {
    si_out[, is_fixed := FALSE]
  }
  
  anno_cols <- intersect(c("gene1", "all_genes", "chr", "start", "end", "gene1_anno"), names(si_out))
  setcolorder(si_out, c("run", "range_id", "event_type", "si", anno_cols))
  si_out[]
}

# Helper function to plot each round
plot_round <- function(round_output) {
  curr_chr <- unique(round_output$chr)
  stopifnot(length(curr_chr) == 1)
  if (!curr_chr %like% "^chr") curr_chr <- paste0("chr", curr_chr)
  
  chr_info <- chr_coordinates[curr_chr, on = "chr_label"][1]
  
  fp <- round_output[
    ,
    .(
      mp          = floor((start + end) / 2)[1],
      gene        = gene1[1],
      si          = median(si),
      is_fixed    = is_fixed[1],
      cancer_anno = gene1_anno[1]
    ),
    by = c("range_id", "event_type")
  ]
  
  fp[cancer_anno %in% c("oncogene", "TSG"),
     gene_label := c(gene, rep(NA, .N - 1)),
     by = c("gene", "event_type")]
  
  palette <- RColorBrewer::brewer.pal(3, "Spectral")
  fp[, inner_fill := fcase(
    cancer_anno == "oncogene", palette[1],
    cancer_anno == "other",    palette[2],
    cancer_anno == "TSG",      palette[3],
    default = "gray"
  )]
  
  fp[is.na(cancer_anno) | cancer_anno == "noncancer", cancer_anno := "noncancer or\nintergenic"]
  palette[4] <- "gray"
  
  fill_legend <- unique(fp[, .(cancer_anno, inner_fill)])
  fill_legend[, legend_order := match(inner_fill, palette)]
  fill_legend <- fill_legend[order(legend_order)]
  
  for_fixed_labels <- fp[is_fixed == TRUE]
  fp[is_fixed == TRUE, gene_label := NA]
  for_fixed_labels[is.na(gene_label), gene_label := gene]
  
  ggplot(fp, aes(x = mp, y = si)) +
    geom_rect(
      data = chr_info,
      aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
      fill = "gray95",
      inherit.aes = FALSE
    ) +
    geom_point(aes(color = inner_fill), size = 2.5) +
    xlab(paste0(curr_chr, " position")) +
    ylab("Scaled selection for copy change") +
    scale_color_identity(
      name   = "Cancer annotation",
      breaks = fill_legend$inner_fill,
      labels = fill_legend$cancer_anno,
      guide  = guide_legend()
    ) +
    facet_wrap(~event_type, ncol = 1, scales = "free") +
    geom_text_repel(
      na.rm              = TRUE,
      aes(label          = gene_label),
      size               = 2,
      min.segment.length = 0,
      point.padding      = 0.5,
      nudge_y            = 0.07,
      segment.colour     = "gray20"
    ) +
    geom_label_repel(
      data               = for_fixed_labels,
      na.rm              = TRUE,
      aes(label          = gene_label),
      size               = 2.5,
      min.segment.length = 0,
      point.padding      = 0.5,
      nudge_y            = 0.07,
      segment.colour     = "gray20",
      color              = "white",
      fill               = "darkmagenta",
      box.padding        = 0.15
    ) +
    scale_x_continuous(labels = scales::label_comma()) +
    theme_light()
}

# This sets our neutral baseline comparison for LRTs
compute_base_loglik <- function(rates, debug = FALSE) {
  etypes <- unique(rates$rates$event_type)
  events_empty <- setNames(replicate(length(etypes), character(0), simplify = FALSE), etypes)
  
  neutral_fit <- cancereffectsizeR:::run_seg_multi(
    events_to_test = events_empty,
    rates          = rates,
    debug          = debug
  )
  
  if ("null_loglik" %in% names(neutral_fit$info)) neutral_fit$info$null_loglik else neutral_fit$info$loglik
}

## Helper functions for bootstrapping
subset_rates_by_samples <- function(rates, samp_vec) {
  dt <- rates$rates
  stopifnot("sample" %in% names(dt))
  
  mult <- as.data.table(table(samp_vec))
  setnames(mult, c("sample", "mult"))
  mult[, mult := as.integer(mult)]
  
  # join to get mult for each sample row, then repeat
  dt_sub <- dt[mult, on = "sample", allow.cartesian = TRUE]
  dt_sub <- dt_sub[rep(seq_len(nrow(dt_sub)), dt_sub$mult)]
  dt_sub[, mult := NULL]
  
  out <- data.table::copy(rates)
  out$rates <- dt_sub
  
  if (!is.null(rates$uncovered) && nrow(rates$uncovered) && "sample" %in% names(rates$uncovered)) {
    uc <- rates$uncovered
    uc_sub <- uc[mult, on = "sample", allow.cartesian = TRUE]
    uc_sub <- uc_sub[rep(seq_len(nrow(uc_sub)), uc_sub$mult)]
    uc_sub[, mult := NULL]
    out$uncovered <- uc_sub
  }
  
  out
}

make_events_list <- function(fixed_dt, cand_row) {
  ev <- rbind(
    fixed_dt[, .(range_id, event_type)],
    cand_row[, .(range_id, event_type)],
    fill = TRUE
  )
  split(ev$range_id, ev$event_type)
}

fit_one_model_loglik <- function(rates, events_to_test, debug = FALSE) {
  fit <- cancereffectsizeR:::run_seg_multi(
    events_to_test = events_to_test,
    rates          = rates,
    debug          = debug
  )
  fit$info$loglik
}

bootstrap_compare_models <- function(
    rates,
    fixed_events,
    cand_A,          # best driver (range_id,event_type)
    cand_B,          # challenger other/noncancer (range_id,event_type)
    B     = 200,
    seed  = 1,
    debug = FALSE,
    cores = 1
) {
  stopifnot(is.data.table(fixed_events), nrow(cand_A) == 1, nrow(cand_B) == 1)
  
  samples <- unique(rates$rates$sample)
  set.seed(seed)
  
  ev_A <- make_events_list(fixed_events, cand_A)
  ev_B <- make_events_list(fixed_events, cand_B)
  
  deltas <- pbapply::pblapply(seq_len(B), function(b) {
    samp_b  <- sample(samples, size = length(samples), replace = TRUE)
    rates_b <- subset_rates_by_samples(rates, samp_b)
    
    llA <- fit_one_model_loglik(rates_b, ev_A, debug = debug)
    llB <- fit_one_model_loglik(rates_b, ev_B, debug = debug)
    llB - llA
  }, cl = cores)
  
  deltas <- unlist(deltas)
  ci <- stats::quantile(deltas, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  p_one_sided <- mean(deltas <= 0, na.rm = TRUE)
  
  list(
    deltas = deltas,
    ci = ci,
    p_one_sided = p_one_sided,
    mean_delta = mean(deltas, na.rm = TRUE)
  )
}

## Iterative selection model with stop determined by LRT
run_iterative_rounds <- function(
    rates,
    output_dir,
    prefix               = "chr",
    cores                = 1,
    base_loglik          = NULL,
    max_rounds           = 30, # for ESCC results we only report the independent results from round 1
    debug                = FALSE,
    
    # The rest of these parameters really only matters for a multi-round application
    alpha_add            = 0.01,  # LRT threshold for adding a new term (nested)
    onco_tsg_annos       = c("oncogene", "TSG"),
    other_annos          = c("other", "noncancer"),
    window_bp            = 2.5e6,   # +/- 2.5 Mb window around best candidate midpoint
    use_bootstrap_override = TRUE,
    B_boot                = 200,
    alpha_boot            = 0.05,
    boot_seed_base        = 123
) {
  fixed_events  <- data.table(range_id = character(), event_type = character())
  round_results <- list()
  
  if (is.null(base_loglik)) {
    message("Computing neutral base_loglik (no CNA effects)...")
    base_loglik <- compute_base_loglik(rates = rates, debug = debug)
    message("  base_loglik = ", signif(base_loglik, 6))
  }
  prev_loglik <- base_loglik
  
  # local helpers
  
  midpoint <- function(dt) floor((dt$start + dt$end) / 2)
  
  choose_within_window <- function(sub, rnd) {
    # sub: candidates within the +/- window, same chr + event_type
    cand_onco  <- sub[gene1_anno %in% onco_tsg_annos]
    cand_other <- sub[gene1_anno %in% other_annos]
    
    best_onco  <- if (nrow(cand_onco))  cand_onco[which.max(loglik)]  else NULL
    best_other <- if (nrow(cand_other)) cand_other[which.max(loglik)] else NULL
    
    # default: prefer driver if present; else best overall
    chosen <- NULL
    if (!is.null(best_onco)) {
      chosen <- best_onco
    } else if (!is.null(best_other)) {
      chosen <- best_other
    } else {
      chosen <- sub[which.max(loglik)]
    }
    
    # if both exist: optionally bootstrap to allow other/noncancer to override driver
    if (!is.null(best_onco) && !is.null(best_other)) {
      ll_gap <- best_other$loglik - best_onco$loglik
      
      if (use_bootstrap_override) {
        if (!is.na(ll_gap) && ll_gap > 0) {
          boot <- bootstrap_compare_models(
            rates        = rates,
            fixed_events = fixed_events,
            cand_A       = best_onco[, .(range_id, event_type)],
            cand_B       = best_other[, .(range_id, event_type)],
            B            = B_boot,
            seed         = boot_seed_base + rnd,
            debug        = debug,
            cores        = cores
          )
          
          ci_low <- unname(boot$ci[[1]])
          p1     <- boot$p_one_sided
          
          if (!is.na(ci_low) && ci_low > 0 && !is.na(p1) && p1 <= alpha_boot) {
            chosen <- best_other
            message(
              "Bootstrap override ACCEPTED: other/noncancer beats driver.",
              " Δloglik(obs)=", signif(ll_gap, 6),
              " boot_CI=[", signif(boot$ci[[1]], 4), ", ", signif(boot$ci[[3]], 4), "]",
              " p_one_sided=", signif(p1, 4)
            )
          } else {
            chosen <- best_onco
            message(
              "Bootstrap override REJECTED: retain driver.",
              " Δloglik(obs)=", signif(ll_gap, 6),
              " boot_CI=[", signif(boot$ci[[1]], 4), ", ", signif(boot$ci[[3]], 4), "]",
              " p_one_sided=", signif(p1, 4)
            )
          }
        } else {
          chosen <- best_onco
          message("Other/noncancer not better by observed loglik; skipping bootstrap and retaining driver.")
        }
      } else {
        # if bootstrap disabled, we simply prefer the best driver by construction
        chosen <- best_onco
      }
    }
    
    chosen
  }
  

  for (rnd in seq_len(max_rounds)) {
    message("=== Round ", rnd, " ===")
    
    r <- if (nrow(fixed_events)) {
      run_round(always_selected = fixed_events, rates = rates, cores = cores, debug = debug)
    } else {
      run_round(rates = rates, cores = cores, debug = debug)
    }
    
    r[, round := rnd]
    round_results[[rnd]] <- r
    
    # Plot
    if (dir.exists(output_dir)) {
      p <- plot_round(r)
      ggsave(
        file.path(output_dir, sprintf("%s_round_%02d.png", prefix, rnd)),
        p, width = 8, height = 6
      )
    }
    
    cand <- r[is_fixed == FALSE & !is.na(loglik)]
    if (!nrow(cand)) {
      message("No unfixed candidates left; stopping.")
      break
    }
    
    # sanity: needed columns for windowing
    needed_cols <- c("chr", "start", "end", "event_type", "gene1", "gene1_anno", "loglik", "range_id")
    missing_cols <- setdiff(needed_cols, names(cand))
    if (length(missing_cols)) {
      stop("Missing columns needed for window-based selection: ", paste(missing_cols, collapse = ", "))
    }
    
    # 1) anchor = best overall candidate by loglik
    anchor <- cand[which.max(loglik)]
    anchor_mp <- midpoint(anchor)
    w_lo <- anchor_mp - window_bp
    w_hi <- anchor_mp + window_bp
    
    # 2) define the “cluster of interest” as +/- window_bp around anchor midpoint
    #    (same chr + same event_type as anchor)
    sub <- cand[
      chr == anchor$chr &
        event_type == anchor$event_type &
        midpoint(cand) >= w_lo &
        midpoint(cand) <= w_hi
    ]
    
    # (if something odd happens and sub is empty, fall back to anchor)
    if (!nrow(sub)) sub <- anchor
    
    message(
      "Anchor (best loglik) at chr=", anchor$chr,
      " event_type=", anchor$event_type,
      " gene=", anchor$gene1,
      " range_id=", anchor$range_id,
      " midpoint=", anchor_mp,
      " window=[", w_lo, ", ", w_hi, "]",
      " candidates_in_window=", nrow(sub),
      " anchor_loglik=", signif(anchor$loglik, 6)
    )
    
    # 3) choose within window using driver preference + bootstrap override
    chosen <- choose_within_window(sub = sub, rnd = rnd)
    
    # Stopping rule: LRT vs current model (df=1)
    ll0 <- prev_loglik
    ll1 <- chosen$loglik
    lrt_stat <- 2 * (ll1 - ll0)
    p_lrt <- stats::pchisq(lrt_stat, df = 1, lower.tail = FALSE)
    
    message(
      "Chosen: gene = ", chosen$gene1,
      " (anno = ", chosen$gene1_anno, ")",
      ", event_type = ", chosen$event_type,
      ", range_id = ", chosen$range_id,
      ", loglik(current) = ", signif(ll0, 6),
      ", loglik(chosen) = ", signif(ll1, 6),
      ", LRT = ", signif(lrt_stat, 6),
      ", p = ", signif(p_lrt, 6)
    )
    
    if (is.na(p_lrt) || p_lrt > alpha_add) {
      message("Stopping: best eligible addition does not pass LRT (p <= ", alpha_add, " required).")
      break
    }
    
    # Accept & fix
    fixed_events <- rbind(
      fixed_events,
      chosen[, .(range_id, event_type)],
      fill = TRUE
    )
    prev_loglik <- ll1
  }
  
  invisible(list(
    fixed_events           = fixed_events,
    round_results          = round_results,
    base_loglik            = base_loglik,
    alpha_add              = alpha_add,
    window_bp              = window_bp,
    use_bootstrap_override = use_bootstrap_override,
    B_boot                 = B_boot,
    alpha_boot             = alpha_boot,
    boot_seed_base         = boot_seed_base
  ))
}


## Run per chromosome
out_dir <- file.path(proj_output_dir, paste0("iterative_", proj, "_chr", chr_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res <- run_iterative_rounds(
  rates                 = chr_rates,
  output_dir            = out_dir,
  prefix                = paste0(proj, "_chr", chr_id),
  cores                 = 1,
  debug                 = FALSE,
  max_rounds            = 30, # our ESCC results reported only focus on round 1
  alpha_add             = 0.01,
  onco_tsg_annos        = c("oncogene", "TSG"),
  other_annos           = c("other", "noncancer"),
  ## bootstrap settings:
  use_bootstrap_override = TRUE,
  B_boot                 = 10000,
  alpha_boot             = 0.05,
  boot_seed_base          = 123
)

saveRDS(res, file.path(out_dir, paste0(proj, "_chr", chr_id, "_iterative_rounds.rds")))

message("Done. Fixed events (head):")
print(head(res$fixed_events))