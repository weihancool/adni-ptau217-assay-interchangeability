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
if (!exists("sb_decision_long", inherits = TRUE))
  stop("Run 03_core_computation.R manually first.", call. = FALSE)
sb_margin_scales <- c("margin_triage", "margin_confirmatory", "margin_youden")
sb_fit_vc <- function(data, outcome, audit, assay_as_fixed = TRUE,
                      include_assay_fold = FALSE, include_assay_run = FALSE) {
  terms <- c("(1|RID)", "(1|occ_id)", "(1|RID:assay)")
  if (!assay_as_fixed) terms <- c(terms, "(1|assay)")
  if (include_assay_fold) terms <- c(terms, "(1|assay_fold)")
  if (include_assay_run)  terms <- c(terms, "(1|assay_run)")
  fixed <- sb_vc_fixed_terms(data, assay_as_fixed)
  f <- as.formula(paste0(
    outcome, " ~ ", paste(fixed, collapse = " + "), " + ",
    paste(terms, collapse = " + ")))
  m <- sb_checked_lmer(f, data = data, REML = TRUE,
                       control = lmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5)))
  sample_counts <- sb_model_sample_counts(m)
  v <- as.data.frame(VarCorr(m))
  getv <- function(g, fitted) {
    if (!fitted) return(NA_real_)
    z <- v$vcov[v$grp == g]
    if (length(z)) z[1] else NA_real_
  }
  vals <- c(person                 = getv("RID", TRUE),
            occasion_within_person = getv("occ_id", TRUE),
            assay                  = getv("assay", !assay_as_fixed),
            person_by_assay        = getv("RID:assay", TRUE),
            assay_by_fold          = getv("assay_fold", include_assay_fold),
            assay_by_run           = getv("assay_run", include_assay_run),
            residual               = getv("Residual", TRUE))
  total <- sum(vals, na.rm = TRUE)
  msg <- m@optinfo$conv$lme4$messages
  tibble(audit = audit, scale = outcome, component = names(vals),
         participants = unname(sample_counts["participants"]),
         occasions = unname(sample_counts["occasions"]),
         n_observations = unname(sample_counts["n_observations"]),
         variance = unname(vals),
         proportion = unname(vals) / total,
         assay_as_fixed = assay_as_fixed,
         singular = isSingular(m),
         converged = sb_lmer_converged(m),
         convergence_message = if (length(msg)) paste(msg, collapse = "; ") else NA_character_,
         person_assay_to_occasion_ratio =
           sb_safe_ratio(vals[["person_by_assay"]], vals[["occasion_within_person"]]),
         ratio_definition = "person_by_assay / shared_occasion_within_person")
}
sb_vc_table3_main <- map_dfr(
  sb_margin_scales,
  ~ sb_fit_vc(sb_decision_long, .x, "crossfit_main")
)
sb_vc_fold <- map_dfr(sb_margin_scales,
                      ~ sb_fit_vc(sb_decision_long, .x, "crossfit_plus_assay_by_fold",
                                  include_assay_fold = TRUE))
sb_run_dat <- sb_decision_long %>% group_by(occ_id) %>%
  filter(all(!is.na(run_no))) %>% ungroup() %>%
  mutate(assay_run = interaction(assay, run_no, drop = TRUE))
sb_vc_run_baseline <- tibble()
sb_vc_run <- if (n_distinct(sb_run_dat$RID) >= 20 &&
                 n_distinct(sb_run_dat$assay_run) > length(study_b_assays)) {
  sb_vc_run_baseline <- map_dfr(sb_margin_scales,
                                ~ sb_fit_vc(sb_run_dat, .x, "crossfit_run_subset_unadjusted"))
  map_dfr(sb_margin_scales,
          ~ sb_fit_vc(sb_run_dat, .x, "crossfit_plus_assay_by_run",
                      include_assay_run = TRUE))
} else {
  message("Insufficient complete collection occasions/RUN_NO levels; assay-by-run audit not assessed.")
  tibble()
}
sb_fixed_reference <- sb_decision_long %>%
  filter(amyloid_positive %in% c(0, 1),
         is.finite(log_value),
         is.finite(abs_amyloid_gap_days),
         !is.na(plasma_date)) %>%
  mutate(assay = as.character(assay)) %>%
  arrange(RID, assay, abs_amyloid_gap_days, plasma_date) %>%
  distinct(RID, assay, .keep_all = TRUE)
sb_fixed_thresholds <- map_dfr(study_b_assays, function(a) {
  d <- filter(sb_fixed_reference, assay == a)
  mu <- mean(d$log_value); sig <- sd(d$log_value)
  z <- (d$log_value - mu) / sig
  th <- sb_pick_thresholds(d$amyloid_positive, z)
  tibble(assay = a, fixed_log_mean = mu, fixed_log_sd = sig,
         fixed_triage_threshold_z = th[["triage"]],
         fixed_confirmatory_threshold_z = th[["confirmatory"]],
         fixed_youden_threshold_z = th[["youden"]],
         reference_participants = n_distinct(d$RID),
         amyloid_positive_n = sum(d$amyloid_positive == 1))
})
sb_write(sb_fixed_thresholds, "Table_S4_fixed_threshold_audit.csv")
sb_decision_fixed <- sb_decision_long %>%
  mutate(assay_chr = as.character(assay)) %>%
  left_join(sb_fixed_thresholds, by = c("assay_chr" = "assay")) %>%
  mutate(fixed_global_z = (log_value - fixed_log_mean) / fixed_log_sd,
         fixed_margin_triage = fixed_global_z - fixed_triage_threshold_z,
         fixed_margin_confirmatory = fixed_global_z - fixed_confirmatory_threshold_z,
         fixed_margin_youden = fixed_global_z - fixed_youden_threshold_z)
sb_vc_fixed <- map_dfr(
  c("fixed_margin_triage", "fixed_margin_confirmatory", "fixed_margin_youden"),
  ~ sb_fit_vc(sb_decision_fixed, .x, "full_sample_fixed_threshold"))
sb_first_available_pet <- sb_decision_long %>%
  arrange(RID, plasma_date, sample_id, assay) %>%
  group_by(RID) %>%
  summarise(
    first_available_amyloid_positive = {
      z <- amyloid_positive[!is.na(amyloid_positive)]
      if (length(z)) first(z) else NA_integer_
    },
    .groups = "drop"
  ) %>%
  mutate(
    baseline_amyloid_group = factor(
      case_when(
        first_available_amyloid_positive == 1L ~ "positive",
        first_available_amyloid_positive == 0L ~ "negative",
        TRUE                                   ~ "unknown"
      ),
      levels = c("negative", "positive", "unknown")
    )
  ) %>%
  select(RID, baseline_amyloid_group)
stopifnot(!anyDuplicated(sb_first_available_pet$RID),
          nrow(sb_first_available_pet) == n_distinct(sb_decision_long$RID))
sb_first_available_dat <- sb_decision_long %>%
  select(-baseline_amyloid_group) %>%
  left_join(sb_first_available_pet, by = "RID", relationship = "many-to-one")
if (anyNA(sb_first_available_dat$baseline_amyloid_group)) {
  stop("Failed to assign first-available PET group to every participant.",
       call. = FALSE)
}
sb_vc_first_available <- sb_fit_vc(sb_first_available_dat,"margin_youden",
                                   "first_available_PET_adjustment")
sb_variance_components_full <- bind_rows(
  sb_vc_table3_main,
  sb_vc_fold,
  sb_vc_run_baseline,
  sb_vc_run,
  sb_vc_fixed,
  sb_vc_first_available
)
sb_variance_audit <- sb_variance_components_full %>%
  filter(scale %in% c("margin_youden", "fixed_margin_youden"),
         component %in% c("person_by_assay", "occasion_within_person",
                          "assay_by_fold", "assay_by_run"),
         !is.na(variance)) %>%
  select(audit, scale, component, variance, proportion, participants, occasions,
         person_assay_to_occasion_ratio, ratio_definition, singular, converged)
sb_write(sb_variance_audit, "Table_S4_person_assay_audit.csv")
cat("\nTable S4 -- person-by-assay three-way audit\n"); print(sb_variance_audit)
sb_write(sb_variance_components_full, "Table_S4_components_full.csv")
