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
if (!sb_have("sb_events", "plasma_meta", "plasma_long", "plasma_wide",
             "amy", "amy_qc_audit", "sb_build_events", "sb_all_events",
             "sb_preprocess_contract_used") ||
    !identical(sb_preprocess_contract_used, sb_make_core_contract())) {
  stop("Run 02_preprocess.R manually first, in this same R session.", call. = FALSE)
}
.sb02_cache <- file.path(study_b_dir, "02_core_results.rds")
.sb02_expected_contract <- sb_make_core_contract()
.sb02_objects <- c("sb_crossfit", "sb_scored", "cohort", "batch_audit",
                   "sb_anchor", "sb_decision_long", "sb_rep_ids", "sb_rep_dat",
                   "sb_events", "sb_core_contract_used")
.sb02_memory_current <- sb_have(.sb02_objects) &&
  identical(sb_core_contract_used, .sb02_expected_contract)
.sb02_needs_recompute <- isTRUE(getOption("study_b_force_rerun", FALSE)) ||
  !.sb02_memory_current || !file.exists(.sb02_cache)
if (.sb02_needs_recompute && file.exists(.sb02_cache) &&
    isTRUE(getOption("study_b_use_cache", FALSE)) &&
    !isTRUE(getOption("study_b_force_rerun"))) {
  message("02_core_computation: loading cached core results from ", .sb02_cache)
  .sb02_cached <- readRDS(.sb02_cache)
  if (!"sb_core_contract_used" %in% names(.sb02_cached) ||
      !identical(.sb02_cached$sb_core_contract_used,
                 .sb02_expected_contract)) {
    stop("02_core_computation: the cached core snapshot was produced from a ",
         "different configuration, code version, or input release. Rerun with ",
         "options(study_b_use_cache = FALSE, study_b_force_rerun = TRUE).",
         call. = FALSE)
  }
  list2env(.sb02_cached, envir = .GlobalEnv)
  if (!sb_have(.sb02_objects)) {
    stop("02_core_computation: cache at ", .sb02_cache, " is incomplete. ",
         "Delete it and rerun.", call. = FALSE)
  }
} else if (.sb02_needs_recompute) {
  sb_core_contract_used <- .sb02_expected_contract
  sb_crossfit <- function(ev) {
    sb_assert_fold_separation(ev)
    stopifnot(!anyDuplicated(ev$sample_id),
              all(vapply(ev[study_b_assays], function(z) all(is.finite(z) & z > 0), logical(1))))
    map_dfr(sort(unique(ev$fold)), function(f) {
      tr <- filter(ev, fold != f)
      te <- filter(ev, fold == f)
      tr_ref <- tr %>% filter(!is.na(amyloid_positive)) %>%
        arrange(RID, abs_amyloid_gap_days, plasma_date) %>%
        distinct(RID, .keep_all = TRUE)
      if (nrow(tr_ref) < 20) {
        stop("Fold ", f, " has only ", nrow(tr_ref),
             " PET-referenced training participants.", call. = FALSE)
      }
      map_dfr(study_b_assays, function(a) {
        log_tr <- log(tr_ref[[a]]); log_te <- log(te[[a]])
        mu <- mean(log_tr); sig <- sd(log_tr)
        if (!is.finite(sig) || sig <= 0) stop("Zero training SD for ", a, call. = FALSE)
        ztr <- (log_tr - mu) / sig
        zte <- (log_te - mu) / sig
        th <- sb_pick_thresholds(tr_ref$amyloid_positive, ztr)
        cal_df <- data.frame(amyloid_positive = tr_ref$amyloid_positive, ztr = ztr)
        cal <- glm(amyloid_positive ~ ztr, data = cal_df, family = binomial())
        if (!isTRUE(cal$converged) || any(!is.finite(coef(cal))) || coef(cal)[2] <= 0) {
          stop("Non-converged or non-monotone calibration model for ", a,
               " in fold ", f, call. = FALSE)
        }
        prob <- as.numeric(stats::predict.glm(cal, newdata = data.frame(ztr = zte),
                                              type = "response"))
        tibble(
          RID = te$RID, sample_id = te$sample_id, plasma_date = te$plasma_date,
          phase = te$phase, occasion = te$occasion,
          years_from_baseline = te$years_from_baseline,
          baseline_amyloid_positive = te$baseline_amyloid_positive,
          baseline_amyloid_group = te$baseline_amyloid_group,
          diagnosis = te$diagnosis,
          primary_anchor = te$primary_anchor,
          amyloid_positive = te$amyloid_positive,
          amyloid_centiloid = te$amyloid_centiloid,
          amyloid_gap_days = te$amyloid_gap_days,
          abs_amyloid_gap_days = te$abs_amyloid_gap_days,
          fold = f, assay = a,
          train_log_mean = mu, train_log_sd = sig,
          training_reference_participants = nrow(tr_ref),
          training_reference_positive_n = sum(tr_ref$amyloid_positive == 1L),
          training_reference_negative_n = sum(tr_ref$amyloid_positive == 0L),
          calibration_boundary = isTRUE(cal$boundary),
          calibration_intercept = unname(coef(cal)[1]),
          calibration_slope = unname(coef(cal)[2]),
          log_value = log_te, z_value = zte, calibrated_probability = prob,
          triage_threshold_z = th[["triage"]],
          confirmatory_threshold_z = th[["confirmatory"]],
          youden_threshold_z = th[["youden"]],
          triage_positive = as.integer(zte >= th[["triage"]]),
          confirmatory_positive = as.integer(zte >= th[["confirmatory"]]),
          youden_positive = as.integer(zte >= th[["youden"]]),
          three_class = factor(case_when(
            zte <  th[["triage"]]       ~ "negative",
            zte >= th[["confirmatory"]] ~ "positive",
            TRUE                        ~ "intermediate"),
            levels = c("negative", "intermediate", "positive")),
          margin_triage = zte - th[["triage"]],
          margin_confirmatory = zte - th[["confirmatory"]],
          margin_youden = zte - th[["youden"]]
        )
      })
    })
  }
  sb_scored <- sb_crossfit(sb_events) %>%
    left_join(plasma_meta %>% select(sample_id, assay, run_no, sample_uid, assay_cv),
              by = c("sample_id", "assay")) %>%
    arrange(RID, plasma_date, assay)
  stopifnot(!anyNA(sb_scored$margin_youden),
            nrow(sb_scored) == nrow(sb_events) * length(study_b_assays),
            !anyDuplicated(sb_scored[, c("sample_id", "assay")]))
  stopifnot(all(table(sb_scored$sample_id) == length(study_b_assays)))
  if (!exists("study_b_has_run_field", inherits = TRUE)) {
    study_b_has_run_field <- any(!is.na(sb_scored$run_no))
  }
  if (!exists("study_b_has_cv_field", inherits = TRUE)) {
    study_b_has_cv_field <- any(!is.na(sb_scored$assay_cv))
  }
  cohort <- tibble(
    source_collection_occasions = plasma_long %>%
      distinct(RID, plasma_date, phase) %>% nrow(),
    complete_four_assay_events = nrow(plasma_wide),
    excluded_incomplete_four_assay_events =
      (plasma_long %>% distinct(RID, plasma_date, phase) %>% nrow()) - nrow(plasma_wide),
    complete_four_assay_participants = n_distinct(plasma_wide$RID),
    amyloid_matched_events =
      n_distinct(sb_scored$sample_id[!is.na(sb_scored$amyloid_positive)]),
    amyloid_matched_participants =
      n_distinct(sb_scored$RID[!is.na(sb_scored$amyloid_positive)]),
    primary_anchor_participants = n_distinct(sb_scored$RID[sb_scored$primary_anchor]),
    anchor_events_with_pet = n_distinct(sb_scored$sample_id[
      sb_scored$primary_anchor & !is.na(sb_scored$amyloid_positive)]),
    anchor_requires_pet = study_b_anchor_requires_pet,
    repeated_participants = sb_scored %>% distinct(RID, sample_id) %>% count(RID) %>%
      summarise(n = sum(n >= 2)) %>% pull(n),
    exact_centiloid_cutoff_scans = sum(amy$amyloid_centiloid ==
                                           study_b_centiloid_cutoff, na.rm = TRUE),
    centiloid_cutoff = study_b_centiloid_cutoff,
    centiloid_positive_operator = ">",
    amyloid_window_days = study_b_amyloid_window,
    pet_qc_allowed = paste(study_b_pet_qc_allowed, collapse = ",")
  )
  .safe_median <- function(x) if (length(x)) stats::median(x) else NA_real_
  batch_audit <- sb_scored %>%
    group_by(assay) %>%
    summarise(
      observations = n(),
      run_no_coverage = mean(!is.na(run_no)),
      sample_id_coverage = mean(!is.na(sample_uid)),
      distinct_runs = n_distinct(run_no[!is.na(run_no)]),
      median_run_size = .safe_median(as.numeric(table(run_no[!is.na(run_no)]))),
      participants_spanning_multiple_runs = sum(tapply(
        run_no, RID, function(z) n_distinct(z[!is.na(z)]) > 1), na.rm = TRUE),
      cv_coverage = mean(!is.na(assay_cv)),
      median_cv = if (any(!is.na(assay_cv))) median(assay_cv, na.rm = TRUE) else NA_real_,
      .groups = "drop") %>%
    mutate(run_field_available = study_b_has_run_field,
           cv_field_available = study_b_has_cv_field)
  if (!study_b_has_run_field) {
    warning("RUN_NO is absent from the FNIH input. The assay-by-run batch audit ",
            "(Table S4) cannot be fitted and must be reported as 'batch identifiers ",
            "not available in the released data', not as a null result.",
            call. = FALSE)
  }
  sb_anchor <- filter(sb_scored, primary_anchor)
  stopifnot(!anyDuplicated(sb_anchor[, c("RID", "assay")]),
            n_distinct(sb_anchor$sample_id) == n_distinct(sb_anchor$RID))
  sb_decision_long <- sb_scored %>%
    mutate(occ_id = interaction(RID, sample_id, drop = TRUE),
           phase = factor(phase), assay = factor(assay, levels = study_b_assays),
           fold = factor(fold),
           assay_fold = interaction(assay, fold, drop = TRUE),
           assay_run = if (study_b_has_run_field)
             interaction(assay, run_no, drop = TRUE) else factor(NA_character_),
           RID_num = RID,
           RID = factor(RID))
  sb_rep_ids <- sb_scored %>% distinct(RID, sample_id) %>% count(RID) %>%
    filter(n >= 2) %>% pull(RID)
  sb_rep_dat <- sb_scored %>% filter(RID %in% sb_rep_ids) %>%
    group_by(RID, sample_id) %>%
    mutate(occasion_mean_margin = mean(margin_youden),
           assay_deviation = margin_youden - occasion_mean_margin,
           assay_spread = sd(margin_youden)) %>%
    ungroup()
  if (nrow(sb_rep_dat) > 0) {
    .dev_fixed <- sb_deviation_fit(sb_rep_dat)
    sb_rep_dat$assay_deviation_adjusted <- as.numeric(residuals(.dev_fixed))
  } else {
    sb_rep_dat$assay_deviation_adjusted <- numeric(0)
    warning("No participant has repeated collection occasions; Figure S6 / Table S9 ",
            "cannot be produced.", call. = FALSE)
  }
  sb_write(sb_scored,   "02_crossfit_decision_scores.csv")
  sb_write(sb_scored %>% distinct(assay, fold, train_log_mean, train_log_sd,
             training_reference_participants, training_reference_positive_n,
             training_reference_negative_n, calibration_boundary,
             calibration_intercept, calibration_slope, triage_threshold_z,
             confirmatory_threshold_z, youden_threshold_z),
           "02_training_fold_audit.csv")
  sb_write(cohort,      "02_cohort_manifest.csv")
  sb_write(batch_audit, "02_batch_field_audit.csv")
  saveRDS(list(sb_crossfit = sb_crossfit, sb_scored = sb_scored, cohort = cohort,
               batch_audit = batch_audit, sb_anchor = sb_anchor,
               sb_decision_long = sb_decision_long, sb_rep_ids = sb_rep_ids,
               sb_rep_dat = sb_rep_dat, sb_events = sb_events,
               sb_core_contract_used = sb_core_contract_used,
               study_b_has_run_field = study_b_has_run_field,
               study_b_has_cv_field = study_b_has_cv_field),
          .sb02_cache)
  cat("\nStudy B cohort\n"); print(cohort)
  cat("\nBatch structure\n"); print(batch_audit)
  cat("\n02_core_computation complete. sb_scored rows:", nrow(sb_scored), "\n")
}
