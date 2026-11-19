if (!exists("sb_theme_classic",mode = "function"))
  stop("Run the updated 01_config.R first in this same R session.",call. = FALSE)
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
if (!exists("sb_begin_reporting_step", mode = "function"))
  stop("Run revised 01_config.R first.", call. = FALSE)
sb_begin_reporting_step("27_headline_synthesis.R")
.headline_core <- sb_read_core()
.sb_need <- function(file, script) {
  p <- file.path(study_b_dir, file)
  if (!file.exists(p))
    stop("27_headline: required output not found: ", file,
         "\nRun ", script, " manually first.", call. = FALSE)
  sb_read_result(p, stringsAsFactors = FALSE)
}
sb_overall_csv   <- .sb_need("Figure_1A_AARR.csv", "07_fig1_interchangeability.R")
sb_gray_csv      <- .sb_need("Figure_1C_three_class_composition.csv", "07_fig1_interchangeability.R")
sb_ratio_ci_csv  <- .sb_need("Figure_S5_variance_ratio_ci.csv", "09_figS5_variance_decomposition.R")
sb_vc_main_csv   <- .sb_need("Figure_S5_variance_components_main.csv", "09_figS5_variance_decomposition.R")
sb_audit_csv     <- .sb_need("Table_S4_person_assay_audit.csv", "10_tableS4_person_assay_audit.R")
sb_persist_csv   <- .sb_need("Table_S9_supporting_longitudinal_persistence.csv", "11_tableS9_longitudinal_summary.R")
sb_pull <- function(df, audit_name, scale_name, component_name, field = "variance") {
  z <- df[df$audit == audit_name & df$scale == scale_name & df$component == component_name, ]
  if (!nrow(z)) return(NA_real_)
  z[[field]][1]
}
sb_pa       <- sb_pull(sb_vc_main_csv, "crossfit_main", "margin_youden", "person_by_assay")
sb_po       <- sb_pull(sb_vc_main_csv, "crossfit_main", "margin_youden", "occasion_within_person")
sb_res      <- sb_pull(sb_vc_main_csv, "crossfit_main", "margin_youden", "residual")
sb_pa_prop  <- sb_pull(sb_vc_main_csv, "crossfit_main", "margin_youden", "person_by_assay", "proportion")
sb_po_prop  <- sb_pull(sb_vc_main_csv, "crossfit_main", "margin_youden", "occasion_within_person", "proportion")
sb_pa_fold      <- sb_pull(sb_audit_csv, "crossfit_plus_assay_by_fold", "margin_youden", "person_by_assay")
sb_po_fold      <- sb_pull(sb_audit_csv, "crossfit_plus_assay_by_fold", "margin_youden", "occasion_within_person")
sb_pa_fold_prop <- sb_pull(sb_audit_csv, "crossfit_plus_assay_by_fold", "margin_youden", "person_by_assay", "proportion")
sb_fold_prop    <- sb_pull(sb_audit_csv, "crossfit_plus_assay_by_fold", "margin_youden", "assay_by_fold", "proportion")
sb_pa_run      <- sb_pull(sb_audit_csv, "crossfit_plus_assay_by_run", "margin_youden", "person_by_assay")
sb_po_run      <- sb_pull(sb_audit_csv, "crossfit_plus_assay_by_run", "margin_youden", "occasion_within_person")
sb_pa_run_prop <- sb_pull(sb_audit_csv, "crossfit_plus_assay_by_run", "margin_youden", "person_by_assay", "proportion")
sb_run_prop    <- sb_pull(sb_audit_csv, "crossfit_plus_assay_by_run", "margin_youden", "assay_by_run", "proportion")
sb_pa_fixed      <- sb_pull(sb_audit_csv, "full_sample_fixed_threshold", "fixed_margin_youden", "person_by_assay")
sb_po_fixed      <- sb_pull(sb_audit_csv, "full_sample_fixed_threshold", "fixed_margin_youden", "occasion_within_person")
sb_pa_fixed_prop <- sb_pull(sb_audit_csv, "full_sample_fixed_threshold", "fixed_margin_youden", "person_by_assay", "proportion")
sb_aarr <- sb_overall_csv$AARR[sb_overall_csv$rule == "youden"]
sb_pa_run_base <- sb_pull(sb_audit_csv, "crossfit_run_subset_unadjusted", "margin_youden", "person_by_assay")
sb_batch_audit_status <- case_when(
  !is.finite(sb_pa_run) ~ "not assessed: assay-by-run model unavailable",
  !is.finite(sb_pa_run_base) ~ "not interpretable: same-subset baseline unavailable",
  TRUE ~ sprintf("assessed; adjusted/unadjusted person-by-assay variance ratio = %.3f",
                 sb_safe_ratio(sb_pa_run, sb_pa_run_base)))
sb_interpretation <- paste(
  if (is.finite(sb_aarr) && sb_aarr > 0)
    "Individual-level cross-assay discordance was observed after cross-fitted calibration."
  else "Interpret cross-assay discordance using the reported estimate and confidence interval.",
  "Variance ratios are descriptive measurement components, not causal effects;",
  "the denominator of the primary ratio is the shared occasion-within-person component."
)
ratio_shared <- sb_ratio_ci_csv[
  sb_ratio_ci_csv$statistic == "person_assay_to_shared_occasion_ratio", , drop = FALSE]
ratio_conservative <- sb_ratio_ci_csv[
  sb_ratio_ci_csv$statistic == "person_assay_to_occasion_plus_residual_ratio", , drop = FALSE]
sb_headline <- tibble(
  primary_anchor_n = n_distinct(.headline_core$sb_anchor$RID),
  anchor_requires_pet = study_b_anchor_requires_pet,
  youden_AARR = sb_aarr,
  youden_AARR_low = sb_overall_csv$AARR_low[sb_overall_csv$rule == "youden"],
  youden_AARR_high = sb_overall_csv$AARR_high[sb_overall_csv$rule == "youden"],
  three_class_disagreement = 1 - sb_gray_csv$percent[
    sb_gray_csv$category == "Complete agreement"] / 100,
  three_class_disagreement_low = 1 - sb_gray_csv$high[
    sb_gray_csv$category == "Complete agreement"] / 100,
  three_class_disagreement_high = 1 - sb_gray_csv$low[
    sb_gray_csv$category == "Complete agreement"] / 100,
  extreme_negative_positive = sb_gray_csv$percent[
    sb_gray_csv$category == "Extreme negative-positive"] / 100,
  extreme_low = sb_gray_csv$low[sb_gray_csv$category == "Extreme negative-positive"] / 100,
  extreme_high = sb_gray_csv$high[sb_gray_csv$category == "Extreme negative-positive"] / 100,
  person_by_assay_variance = sb_pa,
  person_by_assay_proportion = sb_pa_prop,
  shared_occasion_variance = sb_po,
  shared_occasion_proportion = sb_po_prop,
  residual_variance = sb_res,
  person_assay_to_shared_occasion_ratio = sb_safe_ratio(sb_pa, sb_po),
  ratio_low = ratio_shared$low[1], ratio_high = ratio_shared$high[1],
  person_assay_to_occasion_plus_residual_ratio = sb_safe_ratio(sb_pa, sb_po + sb_res),
  conservative_ratio_low = ratio_conservative$low[1],
  conservative_ratio_high = ratio_conservative$high[1],
  audit_assay_by_fold_person_assay_proportion = sb_pa_fold_prop,
  audit_assay_by_fold_component_proportion = sb_fold_prop,
  audit_assay_by_fold_ratio = sb_safe_ratio(sb_pa_fold, sb_po_fold),
  audit_assay_by_run_person_assay_proportion = sb_pa_run_prop,
  audit_assay_by_run_component_proportion = sb_run_prop,
  audit_assay_by_run_ratio = sb_safe_ratio(sb_pa_run, sb_po_run),
  fixed_threshold_person_assay_proportion = sb_pa_fixed_prop,
  fixed_threshold_ratio = sb_safe_ratio(sb_pa_fixed, sb_po_fixed),
  repeated_participants = sb_persist_csv$repeated_participants[1],
  assay_spread_spearman = sb_persist_csv$spread_spearman[1],
  assay_deviation_spearman = sb_persist_csv$assay_deviation_spearman[1],
  assay_deviation_spearman_low = sb_persist_csv$assay_deviation_spearman_low[1],
  assay_deviation_spearman_high = sb_persist_csv$assay_deviation_spearman_high[1],
  assay_deviation_ICC = sb_persist_csv$assay_deviation_ICC[1],
  assay_deviation_ICC_low = sb_persist_csv$assay_deviation_ICC_low[1],
  assay_deviation_ICC_high = sb_persist_csv$assay_deviation_ICC_high[1],
  batch_audit_status = sb_batch_audit_status,
  interpretation = sb_interpretation
)
sb_write(sb_headline, "HEADLINE_results.csv")
cat("\nStudy B headline\n"); print(as.data.frame(sb_headline))
sb_end_reporting_step()
