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
  stop("Cannot locate CODE/", file, ". Open STUDY_B_STANDALONE.Rproj, or set ",
       "Sys.setenv(STUDY_B_ROOT = \"...\").", call. = FALSE)
}
if (!exists("study_b_code_dir", inherits = TRUE))
  stop("Run 01_config.R manually first, in this same R session.", call. = FALSE)
.sb01_cache <- file.path(study_b_dir, "01_preprocessed.rds")
.sb01_expected_contract <- sb_make_core_contract()
.sb01_objects <- c("amy", "amy_pre", "amy_qc_audit", "amy_duplicate_audit",
                   "plasma_duplicate_audit", "plasma_long", "plasma_wide",
                   "plasma_meta", "sb_all_events", "sb_build_events", "sb_events",
                   "sb_preprocess_contract_used", "study_b_has_run_field",
                   "study_b_has_cv_field", "amy_prefilter_audit")
.sb01_memory_current <- sb_have(.sb01_objects) &&
  identical(sb_preprocess_contract_used, .sb01_expected_contract)
.sb01_needs_recompute <- isTRUE(getOption("study_b_force_rerun", FALSE)) ||
  !.sb01_memory_current || !file.exists(.sb01_cache)
if (.sb01_needs_recompute && file.exists(.sb01_cache) &&
    isTRUE(getOption("study_b_use_cache", FALSE)) &&
    !isTRUE(getOption("study_b_force_rerun"))) {
  message("01_preprocess: loading cached preprocessing results from ", .sb01_cache)
  .sb01_cached <- readRDS(.sb01_cache)
  if (!"sb_preprocess_contract_used" %in% names(.sb01_cached) ||
      !identical(.sb01_cached$sb_preprocess_contract_used,
                 .sb01_expected_contract)) {
    stop("01_preprocess: the cached preprocessing snapshot was produced from a ",
         "different configuration, code version, or input release. Rerun with ",
         "options(study_b_use_cache = FALSE, study_b_force_rerun = TRUE).",
         call. = FALSE)
  }
  list2env(.sb01_cached, envir = .GlobalEnv)
  if (!sb_have(.sb01_objects)) {
    stop("01_preprocess: cache at ", .sb01_cache, " is incomplete (missing: ",
         paste(setdiff(.sb01_objects,
                       ls(envir = .GlobalEnv)), collapse = ", "),
         "). Delete it and rerun.", call. = FALSE)
  }
} else if (.sb01_needs_recompute) {
  sb_preprocess_contract_used <- .sb01_expected_contract
  fnih_raw <- sb_load_rda(file.path(study_b_data_dir,
                                    "FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.rda"))
  amy_raw  <- sb_load_rda(file.path(study_b_data_dir, "UCBERKELEY_AMY_6MM.rda"))
  required_fnih_columns <- c("RID", "EXAMDATE", "ORIGPROT", "PLASMA_BIOMARKER", "TESTVALUE")
  required_amy_columns  <- c("RID", "VISCODE", "VISCODE2", "SCANDATE",
                             "CENTILOIDS", "SUMMARY_SUVR", "AMYLOID_STATUS", "qc_flag")
  if (length(setdiff(required_fnih_columns, names(fnih_raw)))) {
    stop("FNIH input is missing required columns: ",
         paste(setdiff(required_fnih_columns, names(fnih_raw)), collapse = ", "))
  }
  if (length(setdiff(required_amy_columns, names(amy_raw)))) {
    stop("Amyloid input is missing required columns: ",
         paste(setdiff(required_amy_columns, names(amy_raw)), collapse = ", "))
  }
  amy_pre <- amy_raw %>%
    transmute(
      RID = as.integer(as.character(RID)),
      VISKEY = sb_viskey(VISCODE, VISCODE2),
      date = sb_safe_date(SCANDATE),
      amyloid_centiloid = sb_safe_num(CENTILOIDS,   sentinel = FALSE),
      amyloid_suvr      = sb_safe_num(SUMMARY_SUVR, sentinel = FALSE),
      amyloid_status    = as.character(AMYLOID_STATUS),
      qc_flag_raw       = suppressWarnings(as.numeric(as.character(qc_flag))),
      qc_pass           = sb_pet_qc_pass(qc_flag)
    ) %>%
    filter(!is.na(RID), !is.na(date),
           is.finite(amyloid_centiloid) | is.finite(amyloid_suvr)) %>%
    distinct()
  amy_eligible <- amy_pre %>%
    filter(qc_pass, is.finite(amyloid_centiloid))
  amy_duplicate_audit <- amy_eligible %>%
    group_by(RID, date) %>%
    summarise(rows = n(),
              distinct_centiloid = n_distinct(
                amyloid_centiloid[is.finite(amyloid_centiloid)]),
              centiloid_range = ifelse(
                sum(is.finite(amyloid_centiloid)) > 0,
                diff(range(amyloid_centiloid[is.finite(amyloid_centiloid)])),
                NA_real_),
              distinct_qc = n_distinct(qc_flag_raw[is.finite(qc_flag_raw)]),
              .groups = "drop") %>%
    filter(rows > 1)
  sb_write(amy_duplicate_audit, "01_amyloid_duplicate_audit.csv")
  if (any(amy_duplicate_audit$distinct_centiloid > 1)) {
    stop("Conflicting QC-eligible amyloid PET rows exist for the same RID/date ",
         "(", sum(amy_duplicate_audit$distinct_centiloid > 1), " participant-dates). ",
         "Resolve them with tracer/processing metadata; see ",
         "01_amyloid_duplicate_audit.csv.", call. = FALSE)
  }
  amy <- amy_eligible %>%
    arrange(RID, date, desc(qc_flag_raw)) %>%
    distinct(RID, date, .keep_all = TRUE)
  qc_raw_counts <- amy_raw %>%
    transmute(qc_flag_raw = suppressWarnings(as.numeric(as.character(qc_flag)))) %>%
    count(qc_flag_raw, name = "scans_raw")
  qc_post_counts <- amy_pre %>%
    count(qc_flag_raw, name = "scans_after_prefilter")
  amy_qc_audit <- full_join(qc_raw_counts, qc_post_counts, by = "qc_flag_raw") %>%
    mutate(scans_raw = coalesce(scans_raw, 0L),
           scans_after_prefilter = coalesce(scans_after_prefilter, 0L),
           scans = scans_after_prefilter,
           eligible = qc_flag_raw %in% study_b_pet_qc_allowed,
           meaning = sb_qc_meaning(qc_flag_raw)) %>%
    arrange(qc_flag_raw)
  amy_prefilter_audit <- tibble(
    step = c("source rows", "after RID/date/measurement filter",
             "after QC filter (qc_flag in allowed set)",
             "after Centiloid-finite filter", "after RID/date de-duplication"),
    rows = c(nrow(amy_raw), nrow(amy_pre),
             sum(amy_pre$qc_pass),
             nrow(amy_eligible), nrow(amy)))
  plasma_long <- fnih_raw %>%
    filter(PLASMA_BIOMARKER %in% unname(study_b_assay_map)) %>%
    transmute(
      RID = as.integer(as.character(RID)),
      plasma_date = sb_safe_date(EXAMDATE),
      phase = coalesce(na_if(toupper(trimws(as.character(ORIGPROT))), ""), "UNKNOWN"),
      assay = names(study_b_assay_map)[match(PLASMA_BIOMARKER, study_b_assay_map)],
      value = sb_safe_num(TESTVALUE, sentinel = FALSE),
      run_no    = if ("RUN_NO"   %in% names(fnih_raw)) sb_clean_id(RUN_NO)   else NA_character_,
      sample_uid= if ("SAMPLEID" %in% names(fnih_raw)) sb_clean_id(SAMPLEID) else NA_character_,
      assay_cv  = if ("CV"       %in% names(fnih_raw)) sb_safe_num(CV, sentinel = FALSE) else NA_real_
    ) %>%
    filter(!is.na(RID), !is.na(plasma_date), !is.na(assay)) %>%
    distinct()
  plasma_duplicate_audit <- plasma_long %>%
    group_by(RID, plasma_date, phase, assay) %>%
    summarise(rows = n(),
              finite_values = n_distinct(value[is.finite(value)]),
              values = paste(sort(unique(signif(value[is.finite(value)], 12))),
                             collapse = ";"),
              runs = paste(sort(unique(run_no[!is.na(run_no)])), collapse = ";"),
              sample_uids = n_distinct(sample_uid[!is.na(sample_uid)]),
              .groups = "drop") %>%
    filter(rows > 1) %>%
    mutate(conflict = finite_values > 1)
  sb_write(plasma_duplicate_audit, "01_plasma_duplicate_audit.csv")
  if (any(plasma_duplicate_audit$sample_uids > 1L))
    stop("Multiple sample identifiers within a RID/date/ORIGPROT/assay cell. ",
         "Do not merge distinct specimens; review 01_plasma_duplicate_audit.csv.", call. = FALSE)
  phase_conflicts <- plasma_long %>% group_by(RID) %>%
    summarise(n_protocols = n_distinct(phase[phase != "UNKNOWN"]), .groups = "drop") %>%
    filter(n_protocols > 1L)
  if (nrow(phase_conflicts))
    stop("ORIGPROT conflicts within RID; resolve enrollment metadata before defining collections.")
  n_conflict <- sum(plasma_duplicate_audit$conflict)
  if (n_conflict > 0) {
    stop(n_conflict, " RID/date/phase/assay cell(s) carry two or more DIFFERENT ",
         "finite values. Resolve with SAMPLEID/RUN_NO/QC; see ",
         "01_plasma_duplicate_audit.csv (conflict == TRUE).", call. = FALSE)
  }
  n_repeat_run <- nrow(plasma_duplicate_audit)
  if (n_repeat_run > 0) {
    message("01_preprocess: ", n_repeat_run, " sample-assay cell(s) appear more ",
            "than once with identical values (repeat runs). Retained; see ",
            "01_plasma_duplicate_audit.csv.")
  }
  plasma_wide <- plasma_long %>%
    group_by(RID, plasma_date, phase, assay) %>%
    summarise(value = sb_first_finite(value), .groups = "drop") %>%
    pivot_wider(names_from = assay, values_from = value) %>%
    filter(if_all(all_of(study_b_assays), ~ is.finite(.x) & .x > 0)) %>%
    arrange(RID, plasma_date) %>%
    mutate(sample_id = row_number())
  if (anyDuplicated(plasma_wide[c("RID", "plasma_date")]))
    stop("Multiple eligible collections for one RID/date. Resolve ORIGPROT/sample metadata; ",
         "same-day duplicates cannot be treated as longitudinal follow-up.", call. = FALSE)
  plasma_meta <- plasma_long %>%
    filter(is.finite(value), value > 0) %>%
    group_by(RID, plasma_date, phase, assay) %>%
    mutate(.unique_run = if (n_distinct(run_no[!is.na(run_no)]) == 1L)
             sb_first_nonmissing(run_no) else NA_character_) %>%
    arrange(is.na(assay_cv), run_no, sample_uid, assay_cv, .by_group = TRUE) %>%
    slice_head(n = 1L) %>% ungroup() %>%
    transmute(RID, plasma_date, phase, assay, run_no = .unique_run,
              sample_uid, assay_cv) %>%
    inner_join(plasma_wide %>% select(RID, plasma_date, phase, sample_id),
               by = c("RID", "plasma_date", "phase"))
  stopifnot(!anyDuplicated(plasma_meta[, c("sample_id", "assay")]))
  study_b_has_run_field <- any(!is.na(plasma_long$run_no))
  study_b_has_cv_field  <- any(!is.na(plasma_long$assay_cv))
  sb_all_events <- plasma_wide %>%
    select(RID, plasma_date, phase, sample_id, all_of(study_b_assays))
  sb_build_events <- function(window_days = study_b_amyloid_window,
                              centiloid_cutoff = study_b_centiloid_cutoff,
                              anchor_requires_pet = study_b_anchor_requires_pet,
                              pet_data = amy,
                              events = sb_all_events,
                              fold_seed = study_b_seed_folds) {
    stopifnot(nrow(events) > 0L, !anyNA(events$sample_id),
              !anyDuplicated(events$sample_id),
              !anyDuplicated(events[c("RID", "plasma_date")]),
              length(window_days) == 1L, is.finite(window_days), window_days > 0,
              length(centiloid_cutoff) == 1L, is.finite(centiloid_cutoff))
    amy_nearest <- events %>%
      inner_join(pet_data %>% select(RID, amy_date = date, amyloid_centiloid,
                                any_of(c("amyloid_status", "VISCODE", "VISCODE2"))),
                 by = "RID", relationship = "many-to-many") %>%
      mutate(amyloid_gap_days = as.numeric(amy_date - plasma_date),
             abs_amyloid_gap_days = abs(amyloid_gap_days)) %>%
      filter(abs_amyloid_gap_days <= window_days, is.finite(amyloid_centiloid)) %>%
      arrange(sample_id, abs_amyloid_gap_days, amy_date) %>%
      distinct(sample_id, .keep_all = TRUE) %>%
      select(sample_id, amy_date, amyloid_centiloid, amyloid_gap_days,
             abs_amyloid_gap_days, any_of(c("amyloid_status", "VISCODE", "VISCODE2")))
    ev <- events %>%
      left_join(amy_nearest, by = "sample_id") %>%
      mutate(amyloid_positive = if_else(
        is.finite(amyloid_centiloid),
        as.integer(amyloid_centiloid > centiloid_cutoff), NA_integer_)) %>%
      arrange(RID, plasma_date) %>%
      group_by(RID) %>%
      mutate(occasion = row_number(),
             baseline_date = first(plasma_date),
             years_from_baseline = as.numeric(plasma_date - baseline_date) / 365.25,
             baseline_amyloid_positive = first(amyloid_positive),
             first_available_amyloid_positive = {
               z <- amyloid_positive[!is.na(amyloid_positive)]
               if (length(z)) first(z) else NA_integer_
             },
             baseline_amyloid_group = factor(case_when(
               baseline_amyloid_positive == 1L ~ "positive",
               baseline_amyloid_positive == 0L ~ "negative",
               TRUE ~ "unknown"), levels = c("negative", "positive", "unknown"))) %>%
      ungroup()
    ev$diagnosis <- NA_character_
    anchors <- ev %>%
      { if (anchor_requires_pet) filter(., !is.na(amyloid_positive)) else . } %>%
      arrange(RID, plasma_date) %>%
      distinct(RID, .keep_all = TRUE) %>%
      select(RID, anchor_sample_id = sample_id)
    ev <- left_join(ev, anchors, by = "RID") %>%
      mutate(primary_anchor = !is.na(anchor_sample_id) & sample_id == anchor_sample_id)
    fold_strata <- ev %>%
      arrange(RID, plasma_date) %>%
      group_by(RID) %>%
      summarise(stratum = {
        z <- amyloid_positive[!is.na(amyloid_positive)]
        if (length(z)) as.character(first(z)) else "unmatched"
      }, .groups = "drop")
    fold_map <- sb_group_folds(fold_strata$RID, fold_strata$stratum, seed = fold_seed)
    left_join(ev, fold_map, by = "RID")
  }
  sb_events <- sb_build_events()
  sb_write(amy_qc_audit, "01_amyloid_qc_audit.csv")
  sb_write(amy_prefilter_audit, "01_amyloid_prefilter_audit.csv")
  saveRDS(list(amy = amy, amy_pre = amy_pre, amy_eligible = amy_eligible,
               amy_qc_audit = amy_qc_audit,
               amy_prefilter_audit = amy_prefilter_audit,
               amy_duplicate_audit = amy_duplicate_audit,
               plasma_duplicate_audit = plasma_duplicate_audit,
               plasma_long = plasma_long, plasma_wide = plasma_wide,
               plasma_meta = plasma_meta, sb_all_events = sb_all_events,
               study_b_has_run_field = study_b_has_run_field,
               study_b_has_cv_field = study_b_has_cv_field,
               sb_build_events = sb_build_events, sb_events = sb_events,
               sb_preprocess_contract_used = sb_preprocess_contract_used),
          .sb01_cache)
  cat("\nAmyloid QC audit (qc_flag is numeric; 0 = fail, -1 = not assessed)\n")
  print(amy_qc_audit)
  cat("\nAmyloid pre-filter accounting\n"); print(amy_prefilter_audit)
  cat("\nBatch metadata available -- RUN_NO: ", study_b_has_run_field,
      " | CV: ", study_b_has_cv_field, "\n", sep = "")
  cat("\n01_preprocess complete. Events built:", nrow(sb_events),
      "| Participants:", n_distinct(sb_events$RID), "\n")
}
