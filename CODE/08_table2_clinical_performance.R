.sb_find <- function(file) {
  if (exists("study_b_code_dir", inherits = TRUE)) return(file.path(study_b_code_dir, file))
  env <- Sys.getenv("STUDY_B_ROOT", unset = "")
  if (nzchar(env) && file.exists(file.path(env, "CODE", file)))
    return(file.path(env, "CODE", file))
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in seq_len(6)) {
    p <- file.path(d, "CODE", file)
    if (file.exists(p)) return(p)
    if (identical(dirname(d), d)) break
    d <- dirname(d)
  }
  stop("Cannot locate CODE/", file, call. = FALSE)
}
if (!exists("study_b_code_dir", inherits = TRUE))
  stop("Run 01_config.R manually first, in this same R session.", call. = FALSE)
if (!sb_have("sb_anchor"))
  stop("Run 02_preprocess.R and 03_core_computation.R manually first.", call. = FALSE)
sb_clinical_performance <- map_dfr(study_b_assays, function(a) {
  g <- sb_anchor %>% filter(assay == a, !is.na(amyloid_positive))
  if (nrow(g) == 0) {
    stop("No anchor observations with a PET reference for assay ", a, call. = FALSE)
  }
  auc_b <- sb_boot_summary(
    g, function(d) sb_auc_rank(d$amyloid_positive, d$calibrated_probability),
    B = study_b_boot_n, seed = sb_seed(paste0("auc_", a)))
  map_dfr(sb_rules, function(col) {
    pred <- g[[col]]; y <- g$amyloid_positive
    rule_lab <- sub("_positive$", "", col)
    sens_b <- sb_boot_summary(
      g, function(d) mean(d[[col]][d$amyloid_positive == 1] == 1),
      B = study_b_boot_n, seed = sb_seed(paste0("sens_", a, "_", rule_lab)))
    spec_b <- sb_boot_summary(
      g, function(d) mean(d[[col]][d$amyloid_positive == 0] == 0),
      B = study_b_boot_n, seed = sb_seed(paste0("spec_", a, "_", rule_lab)))
    tibble(assay = a, rule = rule_lab, n = nrow(g),
           amyloid_positive_n = sum(y == 1),
           amyloid_negative_n = sum(y == 0),
           crossfit_AUC = sb_auc_rank(y, g$calibrated_probability),
           AUC_low = unname(auc_b["low"]), AUC_high = unname(auc_b["high"]),
           AUC_bootstrap_successful = unname(auc_b["bootstrap_successful"]),
           sensitivity = mean(pred[y == 1] == 1),
           sensitivity_low = unname(sens_b["low"]),
           sensitivity_high = unname(sens_b["high"]),
           sensitivity_bootstrap_successful = unname(sens_b["bootstrap_successful"]),
           specificity = mean(pred[y == 0] == 0),
           specificity_low = unname(spec_b["low"]),
           specificity_high = unname(spec_b["high"]),
           specificity_bootstrap_successful = unname(spec_b["bootstrap_successful"]),
           bootstrap_attempted = study_b_boot_n,
           CI_scope = "Participant bootstrap conditional on frozen out-of-fold predictions",
           guideline_scope = "Numerical benchmark only; not clinical guideline validation")
  })
}) %>%
  mutate(meets_AA_triage_90sens_75spec = if_else(
           rule == "triage", sensitivity >= .90 & specificity >= .75, NA),
         meets_AA_confirmatory_90sens_90spec = if_else(
           rule == "confirmatory", sensitivity >= .90 & specificity >= .90, NA))
sb_write(sb_clinical_performance, "Table_2_crossfit_clinical_performance.csv")
.min_success <- min(c(sb_clinical_performance$sensitivity_bootstrap_successful,
                      sb_clinical_performance$specificity_bootstrap_successful,
                      sb_clinical_performance$AUC_bootstrap_successful))
if (.min_success < 0.9 * study_b_boot_n) {
  warning("At least one Table 2 statistic retained only ", .min_success, " of ",
          study_b_boot_n, " bootstrap replicates (degenerate resamples with no ",
          "reference-positive or no reference-negative participants). The ",
          "coverage of the corresponding intervals may be affected -- report the ",
          "retained counts in the supplement.", call. = FALSE)
}
cat("\nTable 2 -- cross-fitted clinical performance\n")
print(sb_clinical_performance)
if (!any(sb_clinical_performance$meets_AA_triage_90sens_75spec, na.rm = TRUE)) {
  cat("\nNOTE: none of the four assays meets the AA triage criterion",
      "(sensitivity >= 90% AND specificity >= 75%) -- report this explicitly",
      "in the main text, not only in the table.\n")
}
