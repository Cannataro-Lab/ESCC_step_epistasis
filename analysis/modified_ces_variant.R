# Function that is nearly identical to ces_variant() in cancereffectsizeR 
# This modified version saves the likelihood function from MLE for each gene to the environment for later calculation of confidence intervals

modified_ces_variant <- function (cesa = NULL, variants = select_variants(cesa, min_freq = 2), 
          samples = character(), model = "default", run_name = "auto", 
          lik_args = list(), optimizer_args = if (identical(model, 
                                                            "default")) list(method = "L-BFGS-B", lower = 0.001, 
                                                                             upper = 1e+09) else list(), return_fit = FALSE, hold_out_same_gene_samples = "auto", 
          cores = 1, conf = 0.95) 
{
  if (!is.numeric(cores) || length(cores) != 1 || cores - 
      as.integer(cores) != 0 || cores < 1) {
    stop("cores should be 1-length positive integer")
  }
  if (!rlang::is_bool(return_fit)) {
    stop("return_fit should be TRUE/FALSE.")
  }
  if (!is(optimizer_args, "list") || uniqueN(names(optimizer_args)) != 
      length(optimizer_args)) {
    stop("optimizer_args should a named list of arguments to pass.")
  }
  reserved_args = c("minuslogl", "start", "vecpar")
  if (any(reserved_args %in% names(optimizer_args))) {
    msg = paste0("Optimizer arguments start, vecpar, and minuslogl cannot be changed here. ", 
                 "If you are using a custom model, your likelihood function can declare these ", 
                 "values directly (see docs).")
    stop(cancereffectsizeR:::pretty_message(msg, emit = F))
  }
  if (!is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  if (length(hold_out_same_gene_samples) == 1) {
    if (!is.logical(hold_out_same_gene_samples)) {
      if (identical(hold_out_same_gene_samples, "auto")) {
        hold_out_same_gene_samples = is(variants, "data.table")
      }
      else {
        stop("hold_out_same_gene_samples should be TRUE/FALSE or left \"auto\".")
      }
    }
  }
  else {
    stop("hold_out_same_gene_samples should be TRUE/FALSE or left \"auto\".")
  }
  if (!is(run_name, "character") || length(run_name) != 1) {
    stop("run_name should be 1-length character")
  }
  if (run_name %in% names(cesa@selection_results)) {
    stop("The run_name you chose has already been used. Please pick a new one.")
  }
  if (!grepl("^[a-z]", tolower(run_name), perl = T) || grepl("\\s\\s", 
                                                             run_name)) {
    stop("Invalid run name. The name must start with a latter and contain no consecutive spaces.")
  }
  if (run_name == "auto") {
    run_number = length(cesa@selection_results) + 1
    run_name = paste0("selection.", run_number)
    while (run_name %in% names(cesa@selection_results)) {
      run_number = run_number + 1
      run_name = paste0("variant_effects_", run_number)
    }
  }
  if (is(model, "character")) {
    model = tolower(model)
    model[model %in% c("sswm", "default")] = "basic"
    model[model %like% "sswm[-_]sequential"] = "sequential"
    if (length(model) != 1 || !model %in% c("basic", "sequential")) {
      stop("model should specify a built-in selection model (i.e., \"default\") or a custom function factory.")
    }
    else {
      if (model == "basic") {
        lik_factory = sswm_lik
      }
      else if (model == "sequential") {
        lik_factory = sswm_sequential_lik
      }
      else {
        stop("Unrecognized model")
      }
    }
  }
  else if (!is(model, "function")) {
    stop("model should specify a built-in selection model (\"default\") or a custom function factory.")
  }
  else {
    lik_factory = model
  }
  if (!is(lik_args, "list")) {
    stop("lik args should be named list")
  }
  if (length(lik_args) != uniqueN(names(lik_args))) {
    stop("lik_args should be a named list without repeated names.")
  }
  samples = cancereffectsizeR:::select_samples(cesa, samples)
  if (samples[, .N] < cesa@samples[, .N]) {
    num_excluded = cesa@samples[, .N] - samples[, .N]
    cancereffectsizeR:::pretty_message(paste0("Note that ", num_excluded, " samples are being excluded from selection inference."))
  }
  cesa = cancereffectsizeR:::copy_cesa(cesa)
  cesa = cancereffectsizeR:::update_cesa_history(cesa, match.call())
  mutations = cesa@mutations
  if (!is.null(conf)) {
    if (is(model, "function")) {
      if (!rlang::is_scalar_double(conf) || conf != 0.95) {
        warning("conf is ignored when running a custom model.")
      }
      conf = NULL
    }
    else {
      if (!is(conf, "numeric") || length(conf) > 1 || 
          conf <= 0 || conf >= 1) {
        stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", 
             call. = F)
      }
    }
  }
  running_compound = FALSE
  if (is(variants, "data.table")) {
    if (!"variant_id" %in% names(variants)) {
      stop("variants table is missing a variant_id column. Typically, variants is generated using select_variants().")
    }
    nonoverlapping = attr(variants, "nonoverlapping")
    if (is.null(nonoverlapping)) {
      if ("variant_id" %in% names(variants)) {
        cancereffectsizeR:::pretty_message("Taking variants from variant_id column of input table....")
      }
    }
    else if (!identical(nonoverlapping, TRUE)) {
      stop("Input variants table may contain overlapping variants; re-run select_variants() to get a non-overlapping table.")
    }
    variants = select_variants(cesa, variant_ids = variants[, 
                                                            variant_id])
  }
  else if (is(variants, "CompoundVariantSet")) {
    running_compound = TRUE
    if (cesa@advanced$uid != variants@cesa_uid) {
      stop("Input CompoundVariantSet does not appear to derive from the input CESAnalysis.")
    }
    if (cesa@samples[, .N] != variants@cesa_num_samples) {
      stop("The number of samples in the CESAnalysis has changed since the CompoundVariantSet was created. ", 
           "Please re-generate it.")
    }
    compound_variants = variants
    variants = select_variants(cesa, variant_ids = compound_variants@snvs$snv_id, 
                               include_subvariants = TRUE)
    variants = variants[compound_variants@snvs, `:=`(compound_name, 
                                                     compound_name), on = c(variant_id = "snv_id")]
    if (variants[, .N] != compound_variants@snvs[, .N]) {
      stop("Internal error: select_variants() didn't return variant info 1-to-1 with compound input.")
    }
    variants[compound_variants@compounds, `:=`(covered_in, 
                                               shared_cov), on = "compound_name"]
  }
  else {
    stop("variants expected to be a variant table (from select_variants(), usually) or a CompoundVariantSet")
  }
  if (variants[, .N] == 0) {
    stop("There are no variants in the input!")
  }
  aac_ids = variants[variant_type == "aac", variant_id]
  noncoding_snv_ids = variants[variant_type == "snv", variant_id]
  if (length(aac_ids) + length(noncoding_snv_ids) == 0) {
    stop("No variants pass filters, so there are no SIs to calculate.", 
         call. = F)
  }
  maf = cesa@maf[samples$Unique_Patient_Identifier, on = "Unique_Patient_Identifier", 
                 nomatch = NULL]
  tmp = unique(maf[, .(gene = unlist(genes)), by = "Unique_Patient_Identifier"])[, 
                                                                                 .(samples = list(Unique_Patient_Identifier)), by = "gene"]
  tumors_with_variants_by_gene = tmp$samples
  names(tumors_with_variants_by_gene) = tmp$gene
  tumors_with_variants_by_gene = list2env(tumors_with_variants_by_gene)
  snv_aac_of_interest = cesa@mutations$aac_snv_key[aac_ids, 
                                                   on = "aac_id"]
  tmp = maf[snv_aac_of_interest, .(variant_id), on = c(variant_id = "snv_id"), 
            by = "Unique_Patient_Identifier", nomatch = NULL]
  tmp[snv_aac_of_interest, `:=`(aac_id, aac_id), on = c(variant_id = "snv_id")]
  tmp = tmp[, .(samples = list(unique(Unique_Patient_Identifier))), 
            by = "aac_id"]
  samples_by_aac = setNames(tmp$samples, tmp$aac_id)
  setkey(maf, "variant_id")
  tmp = maf[noncoding_snv_ids, variant_id, by = "Unique_Patient_Identifier", 
            nomatch = NULL][, .(samples = list(Unique_Patient_Identifier)), 
                            by = "variant_id"]
  samples_by_snv = tmp$samples
  names(samples_by_snv) = tmp$variant_id
  samples_by_variant = list2env(c(samples_by_aac, samples_by_snv))
  setkey(samples, "covered_regions")
  genome_wide_cov_samples = samples["genome", Unique_Patient_Identifier, 
                                    nomatch = NULL]
  selection_results = NULL
  all_coverage = rbind(cesa@mutations$snv[, .(variant_id = snv_id, 
                                              covered_in)], cesa@mutations$amino_acid_change[, .(variant_id = aac_id, 
                                                                                                 covered_in)])
  variants[all_coverage, `:=`(covered_in, covered_in), on = "variant_id"]
  coverage_groups = unique(variants$covered_in)
  num_coverage_groups = length(coverage_groups)
  setkey(maf, "variant_id")
  setkey(variants, "variant_id")
  selection_results = lapply(1:length(coverage_groups), function(i) {
    coverage_group = coverage_groups[[i]]
    if (length(coverage_group) == 1 && is.na(coverage_group)) {
      curr_variants = variants[which(sapply(variants$covered_in, 
                                            function(x) identical(x, NA_character_)))]
    }
    else {
      curr_variants = variants[which(sapply(variants$covered_in, 
                                            function(x) identical(x, coverage_group)))]
    }
    message(sprintf("Preparing to calculate cancer effects (batch %i of %i)...", 
                    i, num_coverage_groups))
    if (is.null(coverage_group)) {
      covered_samples = genome_wide_cov_samples
    }
    else {
      covered_samples = c(samples[coverage_group, Unique_Patient_Identifier, 
                                  nomatch = NULL], genome_wide_cov_samples)
    }
    variants[curr_variants$variant_id, `:=`(num_covered_and_in_samples, 
                                            length(covered_samples)), on = "variant_id"]
    if (length(covered_samples) == 0) {
      warning("Skipped batch ", i, " because no samples had coverage at the variant sites in the batch.")
      next
    }
    work_size = length(covered_samples) * curr_variants[, 
                                                        .N] * 8
    num_proc_groups = ceiling(work_size/1e+09)
    curr_variants[, `:=`(subgroup, ceiling(num_proc_groups * 
                                             1:.N/.N))]
    if (running_compound) {
      curr_variants[, `:=`(subgroup, rep.int(subgroup[1], 
                                             .N)), by = "compound_name"]
      num_proc_groups = max(curr_variants$subgroup)
    }
    curr_results = lapply(1:num_proc_groups, function(j) {
      if (num_proc_groups > 1) {
        message(sprintf("Working on sub-batch %i of %i...", 
                        j, num_proc_groups))
      }
      curr_subgroup = curr_variants[subgroup == j]
      aac_ids = curr_subgroup[variant_type == "aac", variant_id]
      snv_ids = curr_subgroup[variant_type == "snv", variant_id]
      baseline_rates = baseline_mutation_rates(cesa, aac_ids = aac_ids, 
                                               snv_ids = snv_ids, samples = covered_samples)
      gene_lookup = curr_subgroup[, all_genes]
      names(gene_lookup) = curr_subgroup[, variant_id]
      gene_lookup = list2env(gene_lookup)
      process_variant = function(variant_id) {
        if (running_compound) {
          compound_id = variant_id
          tumors_with_variant = intersect(compound_variants@sample_calls[[compound_id]], 
                                          covered_samples)
          current_snvs = compound_variants@snvs[compound_name == 
                                                  compound_id]
          all_genes = current_snvs[, unique(unlist(genes))]
          variant_id = current_snvs$snv_id
          rates = baseline_rates[, ..variant_id]
          rates = rowSums(rates)
        }
        else {
          tumors_with_variant = samples_by_variant[[variant_id]]
          all_genes = gene_lookup[[variant_id]]
          rates = baseline_rates[, ..variant_id][[1]]
        }
        names(rates) = baseline_rates[, Unique_Patient_Identifier]
        if (hold_out_same_gene_samples) {
          if (length(all_genes) == 1) {
            if (is.na(all_genes)) {
              tumors_with_gene_mutated = tumors_with_variant
            }
            else {
              tumors_with_gene_mutated = tumors_with_variants_by_gene[[all_genes]]
            }
          }
          else {
            tumors_with_gene_mutated = unique(unlist(sapply(all_genes, 
                                                            function(x) tumors_with_variants_by_gene[[x]])))
          }
          tumors_without = setdiff(covered_samples, 
                                   tumors_with_gene_mutated)
        }
        else {
          tumors_without = setdiff(covered_samples, 
                                   tumors_with_variant)
        }
        rates_tumors_with = rates[tumors_with_variant]
        rates_tumors_without = rates[tumors_without]
        lik_args = c(list(rates_tumors_with = rates_tumors_with, 
                          rates_tumors_without = rates_tumors_without), 
                     lik_args)
        fn = do.call(lik_factory, lik_args)
        assign(paste0("lik_fn_", all_genes), fn, envir = .GlobalEnv)
        par_init = formals(fn)[[1]]
        names(par_init) = bbmle::parnames(fn)
        final_optimizer_args = c(list(minuslogl = fn, 
                                      start = par_init, vecpar = T), optimizer_args)
        withCallingHandlers({
          fit = do.call(bbmle::mle2, final_optimizer_args)
        }, warning = function(w) {
          if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
            invokeRestart("muffleWarning")
          }
          if (grepl(x = conditionMessage(w), pattern = "convergence failure")) {
            invokeRestart("muffleWarning")
          }
        })
        selection_intensity = bbmle::coef(fit)
        loglikelihood = as.numeric(bbmle::logLik(fit))
        if (running_compound) {
          variant_id = compound_id
        }
        variant_output = c(list(variant_id = variant_id), 
                           as.list(selection_intensity), list(loglikelihood = loglikelihood))
        if (is.character(model) || is.null(lik_args$sample_index)) {
          if (is.character(model)) {
            if (model == "basic") {
              num_samples_with = length(tumors_with_variant)
              num_samples_total = num_samples_with + 
                length(tumors_without)
              variant_output = c(variant_output, list(included_with_variant = num_samples_with, 
                                                      included_total = num_samples_total))
            }
          }
          if (is(model, "function") && is.null(lik_args$sample_index)) {
            num_samples_with = length(tumors_with_variant)
            num_samples_total = num_samples_with + length(tumors_without)
            variant_output = c(variant_output, list(included_with_variant = num_samples_with, 
                                                    included_total = num_samples_total))
          }
        }
        if (!is.null(conf)) {
          min_value = ifelse(is.null(final_optimizer_args$lower), 
                             -Inf, final_optimizer_args$lower)
          max_value = ifelse(is.null(final_optimizer_args$upper), 
                             Inf, final_optimizer_args$upper)
          variant_output = c(variant_output, univariate_si_conf_ints(fit, 
                                                                     fn, min_value, max_value, conf))
        }
        return(list(summary = variant_output, fit = if (return_fit) fit else NULL))
      }
      if (running_compound) {
        variants_to_run = curr_subgroup[, unique(compound_name)]
      }
      else {
        variants_to_run = curr_subgroup$variant_id
      }
      message("Calculating cancer effects...")
      subgroup_results = pbapply::pblapply(variants_to_run, 
                                           process_variant, cl = cores)
      subgroup_selection = rbindlist(lapply(subgroup_results, 
                                            "[[", 1))
      subgroup_fit = lapply(subgroup_results, "[[", 2)
      return(list(subgroup_selection, subgroup_fit))
    })
    group_selection = rbindlist(lapply(curr_results, "[[", 
                                       1))
    group_fit = lapply(curr_results, "[[", 2)
    return(list(group_selection, group_fit))
  })
  fits = unlist(lapply(selection_results, "[[", 2))
  selection_results = rbindlist(lapply(selection_results, 
                                       "[[", 1))
  if (running_compound) {
    selection_results[, `:=`(variant_type, "compound")]
    setnames(selection_results, "variant_id", "variant_name")
    selection_results = selection_results[compound_variants@compounds, 
                                          on = c(variant_name = "compound_name")]
    setattr(selection_results, "is_compound", TRUE)
    setcolorder(selection_results, c("variant_name", "variant_type"))
  }
  else {
    selection_results[variants, `:=`(c("variant_type", "variant_name", 
                                       "gene", "intergenic"), list(variant_type, variant_name, 
                                                                   gene, intergenic)), on = "variant_id"]
    selection_results[intergenic == T, `:=`(gene, NA)]
    selection_results$intergenic = NULL
    setattr(selection_results, "is_compound", FALSE)
    setcolorder(selection_results, c("variant_name", "variant_type", 
                                     "gene"))
    setcolorder(selection_results, c(setdiff(names(selection_results), 
                                             "variant_id"), "variant_id"))
  }
  if (hold_out_same_gene_samples == FALSE) {
    selection_results[, `:=`(held_out, 0)]
  }
  else {
    if ("included_total" %in% names(selection_results)) {
      if (running_compound) {
        num_eligible_by_comp = sapply(compound_variants$definitions, 
                                      function(x) variants[x, min(num_covered_and_in_samples)], 
                                      USE.NAMES = TRUE)
        selection_results[, `:=`(held_out, num_eligible_by_comp[variant_name] - 
                                   included_total)]
      }
      else {
        selection_results[variants, `:=`(held_out, num_covered_and_in_samples - 
                                           included_total), on = "variant_id"]
      }
    }
  }
  if (all(c("held_out", "included_total") %in% names(selection_results))) {
    selection_results[, `:=`(uncovered, samples[, .N] - 
                               included_total - held_out)]
  }
  if (return_fit) {
    fits = lapply(fits, function(x) {
      x@call.orig = call("[not shown]")
      parent.env(environment(x)) = emptyenv()
      return(x)
    })
    setattr(selection_results, "fit", fits)
  }
  curr_results = list(selection_results)
  names(curr_results) = run_name
  cesa@selection_results = c(cesa@selection_results, curr_results)
  return(cesa)
}
