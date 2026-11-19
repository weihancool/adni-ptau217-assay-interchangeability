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
if (!exists("sb_decision_long", inherits = TRUE))
  stop("Run 03_core_computation.R manually first.", call. = FALSE)
.s4_deliver_cols <- c(
  "assay", "fixed_log_mean", "fixed_log_sd",
  "fixed_triage_threshold_z", "fixed_confirmatory_threshold_z",
  "fixed_youden_threshold_z", "reference_participants", "amyloid_positive_n"
)
.s4_figS3_signature <- c(
  "amyloid_negative_n", "triage_log", "confirmatory_log", "youden_log"
)
.s4_is_figS3_product <- function(x) {
  is.data.frame(x) &&
    all(c(.s4_deliver_cols, .s4_figS3_signature) %in% names(x)) &&
    setequal(as.character(x$assay), as.character(study_b_assays)) &&
    !anyDuplicated(as.character(x$assay))
}
.s4_source <- NA_character_
.s4_canonical <- NULL
.s4_in_session <- get0("sb_fixed_thresholds", ifnotfound = NULL)
if (.s4_is_figS3_product(.s4_in_session)) {
  .s4_canonical <- .s4_in_session
  .s4_source <- "in-session object from 20_figS3"
  message("24_tableS5: reusing the fixed-threshold table computed by ",
          "20_figS3 in this session. No second implementation is run.")
} else {
  if (!is.null(.s4_in_session)) {
    message("24_tableS5: an object named sb_fixed_thresholds exists but lacks ",
            "the Fig S3 signature column(s) (",
            paste(setdiff(.s4_figS3_signature, names(.s4_in_session)),
                  collapse = ", "),
            "); it most likely came from 10_tableS4. It will NOT be reused.")
  }
  stop("Run 20_figS3_threshold_distributions.R manually first in this same ",
       "R session. Script 24 will not run it automatically.", call. = FALSE)
}
sb_tableS5_fixed_thresholds <- .s4_canonical %>%
  select(all_of(.s4_deliver_cols)) %>%
  arrange(match(as.character(assay), as.character(study_b_assays)))
stopifnot(
  nrow(sb_tableS5_fixed_thresholds) == length(study_b_assays),
  !anyDuplicated(sb_tableS5_fixed_thresholds$assay),
  identical(as.character(sb_tableS5_fixed_thresholds$assay),
            as.character(study_b_assays))
)
.s4_numeric_cols <- c(
  "fixed_log_mean", "fixed_log_sd",
  "fixed_triage_threshold_z", "fixed_confirmatory_threshold_z",
  "fixed_youden_threshold_z"
)
.s4_nonfinite <- vapply(
  sb_tableS5_fixed_thresholds[.s4_numeric_cols],
  function(z) sum(!is.finite(z)), integer(1)
)
if (any(.s4_nonfinite > 0L)) {
  stop("24_tableS5: non-finite value(s) in ",
       paste(names(.s4_nonfinite)[.s4_nonfinite > 0L], collapse = ", "),
       ". Do not report Table S5 until 20_figS3 is fixed.", call. = FALSE)
}
if (any(sb_tableS5_fixed_thresholds$fixed_log_sd <= 0)) {
  stop("24_tableS5: non-positive fixed_log_sd for assay(s): ",
       paste(sb_tableS5_fixed_thresholds$assay[
         sb_tableS5_fixed_thresholds$fixed_log_sd <= 0], collapse = ", "),
       call. = FALSE)
}
.s4_too_small <- sb_tableS5_fixed_thresholds$reference_participants < 20L
if (any(.s4_too_small)) {
  stop("24_tableS5: too few PET-referenced participants for fixed-threshold ",
       "estimation: ",
       paste(sprintf("%s (n = %d)",
                     sb_tableS5_fixed_thresholds$assay[.s4_too_small],
                     sb_tableS5_fixed_thresholds$reference_participants[.s4_too_small]),
             collapse = ", "),
       call. = FALSE)
}
.s4_bad_order <- with(sb_tableS5_fixed_thresholds,
                      !(fixed_triage_threshold_z <= fixed_youden_threshold_z &
                          fixed_youden_threshold_z <= fixed_confirmatory_threshold_z))
if (any(.s4_bad_order)) {
  warning("24_tableS5: threshold order is not triage <= Youden <= confirmatory ",
          "for: ",
          paste(sb_tableS5_fixed_thresholds$assay[.s4_bad_order], collapse = ", "),
          ". Values are reported as estimated; state this in the table footnote.",
          call. = FALSE)
}
.s4_t3 <- sb_read_result(
  file.path(study_b_dir, "Table_S4_fixed_threshold_audit.csv"),
  check.names = FALSE
)
.s4_compare_cols <- c(.s4_numeric_cols, "reference_participants", "amyloid_positive_n")
if (!all(c("assay", .s4_compare_cols) %in% names(.s4_t3)) ||
    anyDuplicated(.s4_t3$assay))
  stop("Invalid current Table S4 fixed-threshold audit.", call. = FALSE)
.s4_cmp <- inner_join(
  sb_tableS5_fixed_thresholds %>%
    mutate(assay = as.character(assay)) %>%
    select(assay, all_of(.s4_compare_cols)),
  .s4_t3 %>% mutate(assay = as.character(assay)) %>%
    select(assay, all_of(.s4_compare_cols)),
  by = "assay", suffix = c("_S5", "_S4")
)
.s4_diffs <- unlist(lapply(.s4_compare_cols, function(v) {
  abs(.s4_cmp[[paste0(v, "_S5")]] - .s4_cmp[[paste0(v, "_S4")]])
}))
if (nrow(.s4_cmp) != length(study_b_assays) ||
    any(!is.finite(.s4_diffs)) || any(.s4_diffs > 1e-12))
  stop("Table S4 and Table S5 fixed-threshold constants or counts differ.", call. = FALSE)
message("24_tableS5: current S4/S5 fixed-threshold constants and counts agree.")
sb_write(sb_tableS5_fixed_thresholds,
         "Table_S5_fixed_thresholds.csv")
cat("\nTable S5 -- full-sample fixed thresholds (source: ", .s4_source, ")\n",
    sep = "")
print(sb_tableS5_fixed_thresholds)
