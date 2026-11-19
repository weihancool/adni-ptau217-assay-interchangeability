suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
})
N_BOOT <- 1000L
BOOT_SEED <- 20261015L
args <- commandArgs(trailingOnly = TRUE)
find_results <- function() {
  candidates <- c(
    if (length(args)) args[[1]] else character(),
    Sys.getenv("STUDY_B_RESULTS_DIR", unset = ""),
    if (exists("study_b_dir", inherits = TRUE)) get("study_b_dir", inherits = TRUE) else character(),
    if (nzchar(Sys.getenv("STUDY_B_ROOT", unset = "")))
      file.path(Sys.getenv("STUDY_B_ROOT"), "RESULTS") else character(),
    file.path(getwd(), "RESULTS"),
    file.path(dirname(getwd()), "RESULTS")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  found <- candidates[file.exists(file.path(candidates, "02_core_results.rds"))]
  if (!length(found)) {
    stop("Cannot find RESULTS/02_core_results.rds.\n",
         "Run original 01_config.R -> 02_preprocess.R -> 03_core_computation.R first.\n",
         "Or specify RESULTS: Rscript CODE/33_PET_switching_exploratory.R PROJECT/RESULTS\n",
         "Paths checked: ", paste(candidates, collapse = " | "), call. = FALSE)
  }
  normalizePath(found[[1]], winslash = "/", mustWork = TRUE)
}
results_dir <- find_results()
out_dir <- file.path(results_dir, "PET_EXPLORATORY")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("Reading existing core cache: ", file.path(results_dir, "02_core_results.rds"))
core <- readRDS(file.path(results_dir, "02_core_results.rds"))
if (!is.list(core) || is.null(core$sb_scored))
  stop("Cache does not contain core$sb_scored; check original script 03.", call. = FALSE)
d <- as_tibble(core$sb_scored)
need <- c("RID", "sample_id", "plasma_date", "assay", "amyloid_positive",
          "amyloid_centiloid", "abs_amyloid_gap_days", "youden_positive", "three_class")
missing <- setdiff(need, names(d))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
assays <- c("AlzPath", "C2N_ratio", "Fujirebio", "Janssen")
if (!setequal(as.character(unique(d$assay)), assays))
  stop("Unexpected assay names: ", paste(unique(d$assay), collapse = ", "), call. = FALSE)
if (anyDuplicated(d[c("RID", "sample_id", "assay")]))
  stop("Duplicate RID/sample_id/assay rows; cannot construct valid pairs.", call. = FALSE)
if (!all(table(d$sample_id) == 4L))
  stop("Not all specimens have a complete four-assay panel.", call. = FALSE)
d <- d %>%
  mutate(assay = as.character(assay),
         three_class = as.character(three_class),
         plasma_date = as.Date(plasma_date),
         amyloid_positive = as.integer(amyloid_positive),
         youden_positive = as.integer(youden_positive),
         abs_amyloid_gap_days = as.numeric(abs_amyloid_gap_days))
if (any(!is.na(d$amyloid_positive) & !(d$amyloid_positive %in% c(0L, 1L))))
  stop("PET reference must be 0/1 or NA.", call. = FALSE)
if (anyNA(d$youden_positive) ||
    any(!(d$three_class %in% c("negative", "intermediate", "positive"))))
  stop("Incomplete binary or three-class predictions.", call. = FALSE)
if (any(!is.na(d$amyloid_positive) &
        (is.na(d$abs_amyloid_gap_days) | d$abs_amyloid_gap_days > 180)))
  warning("Some PET-matched rows have missing gap or >180 days. ",
          "Check the frozen PET matching configuration.")
fl <- d %>%
  arrange(RID, assay, plasma_date, sample_id) %>%
  group_by(RID, assay) %>%
  summarise(
    n_plasma = n(),
    first_sample = first(sample_id), last_sample = last(sample_id),
    first_date = first(plasma_date), last_date = last(plasma_date),
    first_pet = first(amyloid_positive), last_pet = last(amyloid_positive),
    first_cl = first(amyloid_centiloid), last_cl = last(amyloid_centiloid),
    first_pet_gap = first(abs_amyloid_gap_days),
    last_pet_gap = last(abs_amyloid_gap_days),
    first_binary = first(youden_positive), last_binary = last(youden_positive),
    first_class = first(three_class), last_class = last(three_class),
    .groups = "drop"
  )
audit <- fl %>% group_by(RID) %>% summarise(
  n_assays = n_distinct(assay),
  n_visits_min = min(n_plasma), n_visits_max = max(n_plasma),
  n_first_id = n_distinct(first_sample), n_last_id = n_distinct(last_sample),
  n_first_date = n_distinct(first_date), n_last_date = n_distinct(last_date),
  n_first_pet = n_distinct(first_pet, na.rm = FALSE),
  n_last_pet = n_distinct(last_pet, na.rm = FALSE),
  n_first_pet_gap = n_distinct(first_pet_gap, na.rm = FALSE),
  n_last_pet_gap = n_distinct(last_pet_gap, na.rm = FALSE),
  .groups = "drop"
)
invalid <- audit %>% filter(
  n_assays != 4L | n_visits_min < 2L | n_visits_min != n_visits_max |
  n_first_id != 1L | n_last_id != 1L |
  n_first_date != 1L | n_last_date != 1L |
  n_first_pet != 1L | n_last_pet != 1L |
  n_first_pet_gap != 1L | n_last_pet_gap != 1L
)
if (nrow(invalid)) {
  write.csv(invalid, file.path(out_dir, "PET_AUDIT_invalid_panel.csv"), row.names = FALSE)
  stop("Inconsistent panels: see PET_AUDIT_invalid_panel.csv", call. = FALSE)
}
persons <- fl %>% filter(assay == assays[[1]]) %>%
  mutate(
    years_between = as.numeric(last_date - first_date) / 365.25,
    first_pet_180 = !is.na(first_pet) & !is.na(first_pet_gap) & first_pet_gap <= 180,
    last_pet_180 = !is.na(last_pet) & !is.na(last_pet_gap) & last_pet_gap <= 180,
    first_pet_90 = first_pet_180 & first_pet_gap <= 90,
    last_pet_90 = last_pet_180 & last_pet_gap <= 90,
    both_pet_180 = first_pet_180 & last_pet_180,
    both_pet_90 = first_pet_90 & last_pet_90,
    pet_transition = case_when(
      !both_pet_180 ~ "not both PET-matched",
      first_pet == 0L & last_pet == 0L ~ "negative -> negative",
      first_pet == 0L & last_pet == 1L ~ "negative -> positive",
      first_pet == 1L & last_pet == 0L ~ "positive -> negative",
      first_pet == 1L & last_pet == 1L ~ "positive -> positive",
      TRUE ~ "unexpected"
    )
  ) %>%
  select(RID, first_sample, last_sample, first_date, last_date,
         n_plasma, years_between, first_pet, last_pet, first_cl, last_cl,
         first_pet_gap, last_pet_gap, first_pet_180, last_pet_180,
         first_pet_90, last_pet_90, both_pet_180, both_pet_90, pet_transition)
coverage <- map_dfr(c(180L, 90L), function(w) {
  first_ok <- if (w == 90L) persons$first_pet_90 else persons$first_pet_180
  last_ok <- if (w == 90L) persons$last_pet_90 else persons$last_pet_180
  both_ok <- first_ok & last_ok
  tibble(
    pet_window_days = w,
    all_participants = nrow(persons),
    first_pet_n = sum(first_ok),
    last_pet_n = sum(last_ok),
    both_pet_n = sum(both_ok),
    first_pet_positive_n = sum(first_ok & persons$first_pet == 1L, na.rm = TRUE),
    first_pet_negative_n = sum(first_ok & persons$first_pet == 0L, na.rm = TRUE),
    last_pet_positive_n = sum(last_ok & persons$last_pet == 1L, na.rm = TRUE),
    last_pet_negative_n = sum(last_ok & persons$last_pet == 0L, na.rm = TRUE),
    both_pet_negative_to_positive_n = sum(both_ok & persons$first_pet == 0L &
                                           persons$last_pet == 1L, na.rm = TRUE),
    both_pet_positive_to_negative_n = sum(both_ok & persons$first_pet == 1L &
                                           persons$last_pet == 0L, na.rm = TRUE),
    both_pet_same_status_n = sum(both_ok & persons$first_pet == persons$last_pet,
                                 na.rm = TRUE),
    median_years_between = median(persons$years_between)
  )
})
pet_by_assay <- map_dfr(c(180L, 90L), function(w) {
  map_dfr(c("first", "last"), function(v) {
    rows <- fl %>%
      filter(!is.na(.data[[paste0(v, "_pet")]]),
             !is.na(.data[[paste0(v, "_pet_gap")]]),
             .data[[paste0(v, "_pet_gap")]] <= w) %>%
      transmute(RID, assay, pet = .data[[paste0(v, "_pet")]],
                binary = .data[[paste0(v, "_binary")]],
                class = .data[[paste0(v, "_class")]])
    if (!nrow(rows)) return(tibble())
    rows %>% group_by(assay) %>% summarise(
      pet_window_days = w, visit = v,
      n = n(), pet_positive_n = sum(pet == 1L),
      pet_negative_n = sum(pet == 0L),
      binary_pet_mismatch_per_100 = 100 * mean(binary != pet),
      false_definitive_per_100_all = 100 * mean(
        (pet == 1L & class == "negative") | (pet == 0L & class == "positive")),
      intermediate_per_100 = 100 * mean(class == "intermediate"),
      .groups = "drop"
    )
  })
})
pairs <- map_dfr(assays, function(from) {
  map_dfr(setdiff(assays, from), function(to) {
    s <- fl %>% filter(assay == from)
    t <- fl %>% filter(assay == to)
    joined <- inner_join(s, t, by = "RID", suffix = c("_A", "_B"))
    if (any(joined$first_sample_A != joined$first_sample_B) ||
        any(joined$last_sample_A != joined$last_sample_B))
      stop("Cross-assay blood specimens do not align.", call. = FALSE)
    joined %>% transmute(
      RID, from_assay = from, to_assay = to,
      first_pet = first_pet_A, last_pet = last_pet_A,
      first_pet_gap = first_pet_gap_A, last_pet_gap = last_pet_gap_A,
      first_binary = first_binary_A,
      source_last_binary = last_binary_A,
      target_last_binary = last_binary_B,
      source_last_class = last_class_A,
      target_last_class = last_class_B
    )
  })
}) %>%
  filter(!is.na(last_pet)) %>%
  mutate(
    source_error = as.integer(source_last_binary != last_pet),
    target_error = as.integer(target_last_binary != last_pet),
    error_rd = target_error - source_error,
    harmful = as.integer(source_error == 0L & target_error == 1L),
    beneficial = as.integer(source_error == 1L & target_error == 0L),
    binary_changed = as.integer(source_last_binary != target_last_binary),
    source_false_definitive = as.integer(
      (last_pet == 1L & source_last_class == "negative") |
        (last_pet == 0L & source_last_class == "positive")),
    target_false_definitive = as.integer(
      (last_pet == 1L & target_last_class == "negative") |
        (last_pet == 0L & target_last_class == "positive")),
    false_definitive_rd = target_false_definitive - source_false_definitive,
    source_intermediate = as.integer(source_last_class == "intermediate"),
    target_intermediate = as.integer(target_last_class == "intermediate"),
    three_class_changed = as.integer(source_last_class != target_last_class),
    negative_positive_switch = as.integer(
      (source_last_class == "negative" & target_last_class == "positive") |
        (source_last_class == "positive" & target_last_class == "negative"))
  )
metrics <- c("source_error", "target_error", "error_rd", "harmful", "beneficial",
             "binary_changed", "source_false_definitive", "target_false_definitive",
             "false_definitive_rd", "source_intermediate", "target_intermediate",
             "three_class_changed", "negative_positive_switch")
boot_stats <- function(x, metric_names = metrics, nboot = N_BOOT, seed = BOOT_SEED) {
  if (!nrow(x)) return(tibble(metric = metric_names, per_100 = NA_real_,
                             ci_low = NA_real_, ci_high = NA_real_))
  per_person <- x %>% group_by(RID) %>%
    summarise(across(all_of(metric_names), ~mean(.x)), .groups = "drop")
  X <- as.matrix(per_person[, metric_names, drop = FALSE])
  storage.mode(X) <- "double"
  est <- colMeans(X)
  if (nboot < 1L) {
    lo <- hi <- rep(NA_real_, length(metric_names))
  } else {
    set.seed(seed)
    n <- nrow(X)
    draws <- replicate(nboot, colMeans(X[sample.int(n, n, replace = TRUE), , drop = FALSE]))
    draws <- matrix(draws, nrow = length(metric_names))
    lo <- apply(draws, 1, stats::quantile, probs = 0.025, names = FALSE)
    hi <- apply(draws, 1, stats::quantile, probs = 0.975, names = FALSE)
  }
  tibble(metric = metric_names, per_100 = as.numeric(100 * est),
         ci_low = as.numeric(100 * lo), ci_high = as.numeric(100 * hi))
}
subsets <- list(
  last_pet_180 = filter(pairs, !is.na(last_pet_gap) & last_pet_gap <= 180),
  last_pet_90  = filter(pairs, !is.na(last_pet_gap) & last_pet_gap <= 90),
  both_pet_180 = filter(pairs, !is.na(first_pet), !is.na(first_pet_gap),
                        !is.na(last_pet_gap), first_pet_gap <= 180, last_pet_gap <= 180),
  both_pet_90 = filter(pairs, !is.na(first_pet), !is.na(first_pet_gap),
                       !is.na(last_pet_gap), first_pet_gap <= 90, last_pet_gap <= 90)
)
pair_table <- imap_dfr(subsets, function(dat, nm) {
  if (!nrow(dat)) return(tibble())
  map_dfr(assays, function(f) map_dfr(setdiff(assays, f), function(t) {
    x <- dat %>% filter(from_assay == f, to_assay == t)
    bind_cols(
      tibble(subset = nm, from_assay = f, to_assay = t,
             n_participants = n_distinct(x$RID),
             pet_positive_n = sum(x$last_pet == 1L),
             pet_negative_n = sum(x$last_pet == 0L)),
      boot_stats(x)
    )
  }))
})
source_summary <- imap_dfr(subsets, function(dat, nm) {
  if (!nrow(dat)) return(tibble())
  map_dfr(assays, function(f) {
    x <- dat %>% filter(from_assay == f)
    bind_cols(tibble(subset = nm, from_assay = f,
                     destination = "mean of other 3 assays",
                     n_participants = n_distinct(x$RID),
                     pet_positive_n = n_distinct(x$RID[x$last_pet == 1L]),
                     pet_negative_n = n_distinct(x$RID[x$last_pet == 0L])),
              boot_stats(x))
  })
})
pet_stratified <- map_dfr(c("last_pet_180", "last_pet_90"), function(nm) {
  dat <- subsets[[nm]]
  if (!nrow(dat)) return(tibble())
  map_dfr(c(0L, 1L), function(status) {
    ds <- dat %>% filter(last_pet == status)
    if (!nrow(ds)) return(tibble())
    map_dfr(assays, function(f) {
      x <- ds %>% filter(from_assay == f)
      bind_cols(tibble(subset = nm, PET_status = if (status == 1L) "positive" else "negative",
                       from_assay = f, n_participants = n_distinct(x$RID)),
                boot_stats(x))
    })
  })
})
write_out <- function(x, filename) {
  p <- file.path(out_dir, filename)
  write.csv(x, p, row.names = FALSE, na = "")
  message("Saved: ", p)
}
write_out(coverage, "PET_01_coverage_summary.csv")
write_out(persons, "PET_02_person_level_coverage.csv")
write_out(pet_by_assay, "PET_03_per_assay_first_vs_last.csv")
write_out(pair_table, "PET_04_last_draw_directed_pairwise.csv")
write_out(source_summary, "PET_05_last_draw_by_source_assay.csv")
write_out(pet_stratified, "PET_06_by_source_and_PET_status.csv")
cat("\n==================== PET COVERAGE ====================\n")
print(coverage, width = Inf)
cat("\n================ PRIMARY DIRECTED SWITCHES ================\n")
print(source_summary %>% filter(subset == "last_pet_180",
                                metric %in% c("harmful", "beneficial", "error_rd", "false_definitive_rd")) %>%
        select(from_assay, n_participants, metric, per_100, ci_low, ci_high),
      n = Inf, width = Inf)
cat("\nFinished: results in ", out_dir, "\n", sep = "")
cat("CAUTION: Error means disagreement with temporally matched amyloid PET,\n",
    "NOT proven clinical misdiagnosis. Longitudinal switching here is hypothetical.\n",
    "The combined net risk across all 12 reciprocal directed pairs is 0 by symmetry;\n",
    "interpret specific source->destination comparisons instead.\n", sep = "")
