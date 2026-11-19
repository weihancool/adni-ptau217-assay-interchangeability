suppressPackageStartupMessages({ library(dplyr); library(purrr) })
if (!exists("sb_begin_reporting_step", mode = "function"))
  stop("Run revised 01_config.R first.", call. = FALSE)
sb_begin_reporting_step("31_final_freeze_checks.R")
root <- Sys.getenv("STUDY_B_ROOT")
res <- if (exists("study_b_dir")) study_b_dir else file.path(root, "RESULTS")
core_current <- sb_read_core(file.path(res, "02_core_results.rds"))
read_result <- function(name) sb_read_result(file.path(res, name), check.names = FALSE)
finite_all <- function(x, expected = NULL) {
  length(x) > 0L && (is.null(expected) || length(x) == expected) &&
    !anyNA(x) && all(is.finite(x))
}
date_audit <- read_result("Table_S9_A3_date_alignment_audit.csv")
zero_audit <- read_result("Table_S9_A3_direction_reversal_zero_change_audit.csv")
key_audit <- read_result("Table_S11_B2_sample_key_audit.csv")
b2 <- read_result("Table_S11_B2_variance_explanation.csv")
a2ci <- read_result("Table_S8_A2_PET_performance_exact_CI.csv")
a3res <- read_result("Table_S9_A3_longitudinal_pseudo_switch.csv")
sens <- read_result("Table_S6_sensitivity_window_cutoff_full.csv")
s7 <- read_result("Figure_S7_remeasurement_vs_switching.csv")
a3_all <- function(m) {
  v <- a3res$per_100[a3res$metric == m & a3res$from_assay == "All" & a3res$to_assay == "All"]
  if (length(v) == 1L) v else NA_real_
}
s7_same <- s7$per_100[s7$contrast == "Same assay over time"]
gap_s7 <- if (length(s7_same) == 1L) abs(a3_all("same_assay_three_class_changed") - s7_same) else NA_real_
gap_excess <- abs(a3_all("category_change_excess_due_to_switch") -
                  (a3_all("switched_three_class_changed") - a3_all("same_assay_three_class_changed")))
required_a3_category_metrics <- c(
  "same_assay_three_class_changed",
  "switched_three_class_changed",
  "category_change_excess_due_to_switch"
)
a3_category_metrics_present <- all(required_a3_category_metrics %in% a3res$metric)
sens_anchor_n_present <- "anchor_pet_available_n" %in% names(sens) &&
  any(sens$status == "ok") &&
  all(is.finite(sens$anchor_pet_available_n[sens$status == "ok"]))
joint <- b2 %>% filter(model == "plus_all_assay_by_modifiers")
boot_rate <- if (nrow(joint) && all(c("bootstrap_successful", "bootstrap_attempted") %in% names(joint)))
  joint$bootstrap_successful / joint$bootstrap_attempted else NA_real_
checks <- tibble(
  check = c("A3 has all 12 directed pairs",
            "A3 first dates align within every pair",
            "A3 last dates align within every pair",
            "A3 first sample_id aligns within every pair",
            "A3 last sample_id aligns within every pair",
            "A3 contains no single-draw participants",
            "B2 sample keys are globally safe",
            "A2 exact-CI table has four assays",
            "A3 includes same-assay and switched three-class controls",
            "Sensitivity table reports anchor PET available N",
            "A3 same-assay change equals Figure S7 same-assay estimate",
            "A3 excess equals switched minus same-assay change",
            "B2 joint bootstrap success rate >= 90%"),
  status = c(ifelse(nrow(date_audit) == 12L, "PASS", "FAIL"),
             ifelse(finite_all(date_audit$pct_first_date_match, 12L) && all(date_audit$pct_first_date_match == 100), "PASS", "FAIL"),
             ifelse(finite_all(date_audit$pct_last_date_match, 12L) && all(date_audit$pct_last_date_match == 100), "PASS", "FAIL"),
             ifelse(finite_all(date_audit$pct_first_sample_match, 12L) && all(date_audit$pct_first_sample_match == 100), "PASS", "FAIL"),
             ifelse(finite_all(date_audit$pct_last_sample_match, 12L) && all(date_audit$pct_last_sample_match == 100), "PASS", "FAIL"),
             ifelse(finite_all(date_audit$n_single_draw, 12L) && all(date_audit$n_single_draw == 0), "PASS", "FAIL"),
             ifelse(nrow(key_audit) == 6L && !anyNA(key_audit$passed) && all(key_audit$passed), "PASS", "FAIL"),
             ifelse(nrow(a2ci) == 4L && setequal(a2ci$assay, study_b_assays) &&
                    all(a2ci$pet_positive_n > 0 & a2ci$pet_negative_n > 0), "PASS", "FAIL"),
             ifelse(a3_category_metrics_present, "PASS", "FAIL"),
             ifelse(sens_anchor_n_present, "PASS", "FAIL"),
             ifelse(is.finite(gap_s7) && gap_s7 < 1e-8, "PASS", "FAIL"),
             ifelse(is.finite(gap_excess) && gap_excess < 1e-8, "PASS", "FAIL"),
             ifelse(is.finite(boot_rate) && boot_rate >= .90, "PASS", "WARN")),
  value = c(nrow(date_audit), min(date_audit$pct_first_date_match),
            min(date_audit$pct_last_date_match), min(date_audit$pct_first_sample_match),
            min(date_audit$pct_last_sample_match), sum(date_audit$n_single_draw),
            sum(key_audit$passed), nrow(a2ci),
            sum(required_a3_category_metrics %in% a3res$metric),
            if ("anchor_pet_available_n" %in% names(sens)) {
              ok_n <- sens$anchor_pet_available_n[sens$status == "ok"]
              if (length(ok_n)) min(ok_n, na.rm = TRUE) else NA_real_
            } else NA_real_,
            gap_s7, gap_excess,
            boot_rate),
  note = c("Hard design check", "Hard same-draw check", "Hard same-draw check",
           "Hard same-sample check", "Hard same-sample check", "Hard repeated-draw check",
           "Hard join-key check", "Reporting completeness",
           "A3 definition completeness", "Sensitivity sample-size transparency",
           "Cross-script identity (value = absolute difference, per 100)",
           "Internal identity (value = absolute difference, per 100)",
           "Exploratory variance-component bootstrap")
)
legacy_diff <- zero_audit %>% filter(from_assay == "All", to_assay == "All") %>%
  summarise(n_definition_discordant = sum(n_definition_discordant, na.rm = TRUE),
            n_delta_same_zero = sum(n_delta_same_zero, na.rm = TRUE),
            n_delta_switch_zero = sum(n_delta_switch_zero, na.rm = TRUE))
sb_save_csv(checks, file.path(res, "FINAL_FREEZE_CHECKS.csv"), row.names = FALSE, na = "")
sb_save_csv(legacy_diff, file.path(res, "FINAL_FREEZE_A3_STRICT_VS_LEGACY.csv"), row.names = FALSE, na = "")
hard_fail <- checks$status == "FAIL"
if (anyNA(checks$status) || any(hard_fail)) stop("Final freeze audit failed: ", paste(checks$check[hard_fail], collapse = "; "), call. = FALSE)
if (any(checks$status == "WARN")) warning(paste(checks$check[checks$status == "WARN"], collapse = "; "))
message("Final freeze hard checks passed.")
sb_assert_fold_separation(core_current$sb_events)
stopifnot(!anyDuplicated(core_current$sb_events[c("RID", "plasma_date")]),
          !anyDuplicated(core_current$sb_scored[c("sample_id", "assay")]))
.aarr <- sb_read_result("Figure_1A_AARR.csv", check.names = FALSE)
if (any(!is.finite(.aarr$AARR_low)) || any(!is.finite(.aarr$AARR_high)))
  stop("Primary AARR intervals missing/nonfinite; freeze is incomplete.")
if (any(.aarr$pipeline_bootstrap_successful < .90 * .aarr$pipeline_bootstrap_attempted))
  warning("Primary AARR bootstrap success below 90%; inspect failures before submission.")
.model_rows <- get0("sb_model_audit", ifnotfound = NULL)
if (is.environment(.model_rows) && length(.model_rows$rows))
  sb_write(bind_rows(.model_rows$rows), "FINAL_REPORTING_MODEL_DIAGNOSTICS.csv")
writeLines(capture.output(sessionInfo()), file.path(res, "FINAL_SESSION_INFO.txt"))
message("Current-core identity, result hashes, and primary CI checks passed.")
.s4_final <- read_result("Table_S4_variance_components_full.csv")
.s4_counts <- read_result("Table_S4_model_sample_counts_audit.csv")
.s4_keys <- .s4_final %>% distinct(audit, scale)
if (anyDuplicated(.s4_counts[c("audit","scale")]) ||
    !setequal(paste(.s4_keys$audit,.s4_keys$scale),
              paste(.s4_counts$audit,.s4_counts$scale)) ||
    any(!is.finite(as.matrix(.s4_counts[c("participants","occasions","n_observations")]))) ||
    any(as.matrix(.s4_counts[c("participants","occasions","n_observations")]) <= 0) ||
    !all(.s4_counts$variance_match %in% TRUE))
  stop("Verified S4 model metadata incomplete. Run 32_reporting_input_audit.R.", call. = FALSE)
sb_write(.s4_counts,"FINAL_S4_MODEL_SAMPLE_COUNTS.csv")
sb_finalize_display_outputs(res)
message("Supplementary Table S1-S11 section bundles, display map and export checks written.")
sb_end_reporting_step()
