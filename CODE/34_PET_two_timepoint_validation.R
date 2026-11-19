suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
})
N_BOOT <- 1000L
BOOT_SEED <- 20261015L
MIN_N_FOR_CI <- 10L
CL_CUTOFF <- 20
CL_EXCLUSION_MARGIN <- 5
ASSAYS <- c("AlzPath", "C2N_ratio", "Fujirebio", "Janssen")
find_results <- function() {
  args <- commandArgs(trailingOnly = TRUE)
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
  if (!length(found)) stop(
    "Cannot locate 02_core_results.rds. Set STUDY_B_RESULTS_DIR to your RESULTS folder.\nChecked: ",
    paste(candidates, collapse = " | "), call. = FALSE)
  normalizePath(found[[1]], winslash = "/", mustWork = TRUE)
}
results_dir <- find_results()
out_dir <- file.path(results_dir, "PET_LONGITUDINAL_EXPLORATORY")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Reading: ", file.path(results_dir, "02_core_results.rds"))
core <- readRDS(file.path(results_dir, "02_core_results.rds"))
if (!is.list(core) || is.null(core$sb_scored))
  stop("Core cache lacks $sb_scored. Run your original script 03 first.", call. = FALSE)
d <- as_tibble(core$sb_scored)
required <- c("RID", "sample_id", "plasma_date", "assay", "amyloid_positive",
              "amyloid_centiloid", "amyloid_gap_days", "abs_amyloid_gap_days",
              "youden_positive", "three_class")
missing_columns <- setdiff(required, names(d))
if (length(missing_columns))
  stop("Required columns missing: ", paste(missing_columns, collapse = ", "), call. = FALSE)
d <- d %>% mutate(
  assay = as.character(assay), plasma_date = as.Date(plasma_date),
  amyloid_positive = as.integer(amyloid_positive),
  amyloid_centiloid = as.numeric(amyloid_centiloid),
  amyloid_gap_days = as.numeric(amyloid_gap_days),
  abs_amyloid_gap_days = as.numeric(abs_amyloid_gap_days),
  youden_positive = as.integer(youden_positive),
  three_class = as.character(three_class)
)
if (!setequal(unique(d$assay), ASSAYS)) stop("Assay names differ from original study.")
if (anyNA(d$RID) || anyNA(d$sample_id) || anyNA(d$plasma_date))
  stop("Missing RID/sample_id/plasma_date.")
if (anyDuplicated(d[c("RID", "sample_id", "assay")]))
  stop("Duplicate sample/assay records: revisit core data.")
if (any(table(d$sample_id) != length(ASSAYS)))
  stop("Every sample_id must have a complete four-assay panel.")
if (anyNA(d$youden_positive) ||
    any(!d$youden_positive %in% 0:1) ||
    anyNA(d$three_class) ||
    any(!d$three_class %in% c("negative", "intermediate", "positive")))
  stop("Missing or invalid cross-fitted binary / three-level classifications.")
if (any(!is.na(d$amyloid_positive) & !d$amyloid_positive %in% 0:1))
  stop("PET status must be 0, 1, or NA.")
if (any(!is.na(d$amyloid_positive) &
        (!is.finite(d$amyloid_centiloid) | is.na(d$abs_amyloid_gap_days))))
  stop("PET-positive/negative entries must have finite Centiloid and gap.")
if (any(!is.na(d$amyloid_positive) &
        d$amyloid_positive != as.integer(d$amyloid_centiloid > CL_CUTOFF)))
  stop("PET classifications disagree with Centiloid >20; inspect configuration.")
fl <- d %>% arrange(RID, assay, plasma_date, sample_id) %>%
  group_by(RID, assay) %>% summarise(
    n_draws = n(),
    first_sample = first(sample_id), last_sample = last(sample_id),
    first_date = first(plasma_date), last_date = last(plasma_date),
    first_pet = first(amyloid_positive), last_pet = last(amyloid_positive),
    first_cl = first(amyloid_centiloid), last_cl = last(amyloid_centiloid),
    first_gap = first(abs_amyloid_gap_days), last_gap = last(abs_amyloid_gap_days),
    first_signed_gap = first(amyloid_gap_days),
    last_signed_gap = last(amyloid_gap_days),
    first_binary = first(youden_positive), last_binary = last(youden_positive),
    first_class = first(three_class), last_class = last(three_class),
    .groups = "drop"
  )
audit <- fl %>% group_by(RID) %>% summarise(
  n_assays = n_distinct(assay), n_min = min(n_draws), n_max = max(n_draws),
  first_samples = n_distinct(first_sample, na.rm = FALSE),
  last_samples = n_distinct(last_sample, na.rm = FALSE),
  first_dates = n_distinct(first_date, na.rm = FALSE),
  last_dates = n_distinct(last_date, na.rm = FALSE),
  first_pets = n_distinct(first_pet, na.rm = FALSE),
  last_pets = n_distinct(last_pet, na.rm = FALSE),
  first_cls = n_distinct(first_cl, na.rm = FALSE),
  last_cls = n_distinct(last_cl, na.rm = FALSE),
  first_gaps = n_distinct(first_gap, na.rm = FALSE),
  last_gaps = n_distinct(last_gap, na.rm = FALSE),
  .groups = "drop"
)
bad <- audit %>% filter(
  n_assays != 4L | n_min < 2L | n_min != n_max |
  first_samples != 1L | last_samples != 1L |
  first_dates != 1L | last_dates != 1L |
  first_pets != 1L | last_pets != 1L |
  first_cls != 1L | last_cls != 1L |
  first_gaps != 1L | last_gaps != 1L
)
if (nrow(bad)) {
  write.csv(bad, file.path(out_dir, "34_AUDIT_inconsistent_panel.csv"), row.names = FALSE)
  stop("First/last panels disagree between assays. Inspect 34_AUDIT_inconsistent_panel.csv")
}
persons <- fl %>% filter(assay == ASSAYS[[1]]) %>% transmute(
  RID, first_date, last_date, n_draws,
  years_between = as.numeric(last_date - first_date) / 365.25,
  first_pet, last_pet, first_cl, last_cl, first_gap, last_gap,
  first_pet_date = first_date + first_signed_gap,
  last_pet_date = last_date + last_signed_gap,
  same_pet_date_at_both_ends = !is.na(first_signed_gap) & !is.na(last_signed_gap) &
    first_pet_date == last_pet_date,
  near_cutoff_at_either = !is.na(first_cl) & !is.na(last_cl) &
    (abs(first_cl - CL_CUTOFF) <= CL_EXCLUSION_MARGIN |
       abs(last_cl - CL_CUTOFF) <= CL_EXCLUSION_MARGIN),
  delta_centiloid = last_cl - first_cl,
  pet_transition = case_when(
    is.na(first_pet) | is.na(last_pet) ~ "not_both_pet_matched",
    first_pet == 0L & last_pet == 1L ~ "PET_0_to_1",
    first_pet == 0L & last_pet == 0L ~ "PET_0_to_0",
    first_pet == 1L & last_pet == 1L ~ "PET_1_to_1",
    first_pet == 1L & last_pet == 0L ~ "PET_1_to_0",
    TRUE ~ "unexpected"
  )
)
write_out <- function(x, filename) {
  file <- file.path(out_dir, filename)
  utils::write.csv(x, file = file, row.names = FALSE, na = "")
  message("Saved: ", file)
}
write_out(persons, "34_01_person_pet_transition_audit.csv")
pairs <- map_dfr(ASSAYS, function(from) map_dfr(setdiff(ASSAYS, from), function(to) {
  a <- fl %>% filter(assay == from)
  b <- fl %>% filter(assay == to)
  ab <- inner_join(a, b, by = "RID", suffix = c("_A", "_B"))
  if (nrow(ab) != nrow(persons) ||
      any(ab$first_sample_A != ab$first_sample_B) ||
      any(ab$last_sample_A != ab$last_sample_B))
    stop("Invalid matching between two assays at first/last blood samples.")
  ab %>% transmute(
    RID, from_assay = from, to_assay = to,
    first_date = first_date_A, last_date = last_date_A,
    first_pet = first_pet_A, last_pet = last_pet_A,
    first_cl = first_cl_A, last_cl = last_cl_A,
    first_gap = first_gap_A, last_gap = last_gap_A,
    A_first = first_binary_A, A_last = last_binary_A,
    B_first = first_binary_B, B_last = last_binary_B,
    A_first_class = first_class_A, A_last_class = last_class_A,
    B_first_class = first_class_B, B_last_class = last_class_B
  )
})) %>% mutate(
  pet_transition = case_when(
    first_pet == 0L & last_pet == 1L ~ "PET_0_to_1",
    first_pet == 0L & last_pet == 0L ~ "PET_0_to_0",
    first_pet == 1L & last_pet == 1L ~ "PET_1_to_1",
    first_pet == 1L & last_pet == 0L ~ "PET_1_to_0",
    TRUE ~ "not_both_pet_matched"
  ),
  source_first_pet_consistent = !is.na(first_pet) & A_first == first_pet,
  first_A_error = as.integer(A_first != first_pet),
  first_B_error = as.integer(B_first != first_pet),
  last_A_error = as.integer(A_last != last_pet),
  last_B_error = as.integer(B_last != last_pet),
  first_error_rd = first_B_error - first_A_error,
  last_error_rd = last_B_error - last_A_error,
  delta_error_rd = last_error_rd - first_error_rd,
  source_last_harmful = as.integer(last_A_error == 0L & last_B_error == 1L),
  source_last_beneficial = as.integer(last_A_error == 1L & last_B_error == 0L),
  A_last_positive = as.integer(A_last == 1L),
  B_last_positive = as.integer(B_last == 1L),
  positive_rd = B_last_positive - A_last_positive,
  A_apparent_up = as.integer(A_first == 0L & A_last == 1L),
  B_apparent_up = as.integer(A_first == 0L & B_last == 1L),
  apparent_up_rd = B_apparent_up - A_apparent_up,
  A_apparent_down = as.integer(A_first == 1L & A_last == 0L),
  B_apparent_down = as.integer(A_first == 1L & B_last == 0L),
  apparent_down_rd = B_apparent_down - A_apparent_down,
  A_sequence_pet_consistent = as.integer(A_first == first_pet & A_last == last_pet),
  B_sequence_pet_consistent = as.integer(A_first == first_pet & B_last == last_pet),
  sequence_consistency_rd = B_sequence_pet_consistent - A_sequence_pet_consistent,
  A_false_definitive = as.integer(
    (last_pet == 1L & A_last_class == "negative") |
      (last_pet == 0L & A_last_class == "positive")),
  B_false_definitive = as.integer(
    (last_pet == 1L & B_last_class == "negative") |
      (last_pet == 0L & B_last_class == "positive")),
  false_definitive_rd = B_false_definitive - A_false_definitive,
  last_pair_reclassified = as.integer(A_last != B_last),
  informative_first = as.integer(abs(first_error_rd) == 1L),
  informative_last = as.integer(abs(last_error_rd) == 1L),
  informative_both = as.integer(informative_first == 1L & informative_last == 1L),
  same_favored_when_informative = as.integer(informative_both == 1L &
                                                first_error_rd == last_error_rd)
)
specs <- tibble(
  subset = c("PET90", "PET180", "PET90_no_CL15to25", "PET180_no_CL15to25"),
  window = c(90L, 180L, 90L, 180L),
  exclude_near_cutoff = c(FALSE, FALSE, TRUE, TRUE)
)
pet_filter <- function(x, window, exclude_near_cutoff = FALSE) {
  y <- x %>% filter(
    !is.na(first_pet), !is.na(last_pet),
    !is.na(first_gap), !is.na(last_gap),
    first_gap <= window, last_gap <= window,
    is.finite(first_cl), is.finite(last_cl)
  )
  if (exclude_near_cutoff) y <- y %>% filter(
    abs(first_cl - CL_CUTOFF) > CL_EXCLUSION_MARGIN,
    abs(last_cl - CL_CUTOFF) > CL_EXCLUSION_MARGIN
  )
  y
}
coverage <- map_dfr(seq_len(nrow(specs)), function(i) {
  s <- specs[i, ]
  y <- persons %>% filter(
    !is.na(first_pet), !is.na(last_pet),
    !is.na(first_gap), !is.na(last_gap),
    first_gap <= s$window, last_gap <= s$window,
    is.finite(first_cl), is.finite(last_cl)
  )
  if (s$exclude_near_cutoff) y <- y %>% filter(
    abs(first_cl - CL_CUTOFF) > CL_EXCLUSION_MARGIN,
    abs(last_cl - CL_CUTOFF) > CL_EXCLUSION_MARGIN
  )
  y %>% summarise(
    subset = s$subset, n_both_pet = n(),
    n_pet_0_to_1 = sum(pet_transition == "PET_0_to_1"),
    n_pet_0_to_0 = sum(pet_transition == "PET_0_to_0"),
    n_pet_1_to_1 = sum(pet_transition == "PET_1_to_1"),
    n_pet_1_to_0 = sum(pet_transition == "PET_1_to_0"),
    n_same_pet_date_both = sum(same_pet_date_at_both_ends),
    n_near_cutoff_either = sum(near_cutoff_at_either),
    median_years_between = median(years_between),
    .groups = "drop"
  )
})
write_out(coverage, "34_02_pet_coverage_and_transitions.csv")
METRICS <- c(
  "last_A_error", "last_B_error", "last_error_rd",
  "source_last_harmful", "source_last_beneficial",
  "A_last_positive", "B_last_positive", "positive_rd",
  "A_apparent_up", "B_apparent_up", "apparent_up_rd",
  "A_apparent_down", "B_apparent_down", "apparent_down_rd",
  "A_sequence_pet_consistent", "B_sequence_pet_consistent", "sequence_consistency_rd",
  "A_false_definitive", "B_false_definitive", "false_definitive_rd",
  "last_pair_reclassified"
)
PERSISTENCE_METRICS <- c(
  "first_A_error", "first_B_error", "first_error_rd",
  "last_A_error", "last_B_error", "last_error_rd", "delta_error_rd",
  "informative_first", "informative_last", "informative_both",
  "same_favored_when_informative"
)
boot_metrics <- function(dat, cols, seed = BOOT_SEED) {
  if (!nrow(dat)) return(tibble(
    metric = cols, per_100 = NA_real_, ci_low = NA_real_, ci_high = NA_real_
  ))
  z <- dat %>% group_by(RID) %>% summarise(
    across(all_of(cols), ~ mean(.x, na.rm = FALSE)), .groups = "drop"
  )
  m <- as.matrix(z[, cols, drop = FALSE]); storage.mode(m) <- "double"
  if (anyNA(m)) stop("Unexpected missing metric among PET-matched rows.")
  point <- colMeans(m)
  lower <- upper <- rep(NA_real_, length(cols))
  if (nrow(m) >= MIN_N_FOR_CI && N_BOOT > 0L) {
    set.seed(seed)
    n_id <- nrow(m)
    resample_counts <- stats::rmultinom(n = N_BOOT, size = n_id,
                                        prob = rep(1 / n_id, n_id))
    boot <- (t(m) %*% resample_counts) / n_id
    lower <- apply(boot, 1L, stats::quantile, probs = 0.025, names = FALSE)
    upper <- apply(boot, 1L, stats::quantile, probs = 0.975, names = FALSE)
  }
  tibble(metric = cols, per_100 = 100 * as.numeric(point),
         ci_low = 100 * lower, ci_high = 100 * upper)
}
make_table <- function(dat, cols, by_source = FALSE, seed = BOOT_SEED) {
  if (!nrow(dat)) return(tibble())
  if (by_source) {
    map_dfr(ASSAYS, function(from) {
      x <- dat %>% filter(from_assay == from)
      if (!nrow(x)) return(tibble())
      bind_cols(
        tibble(from_assay = from, to_assay = "mean_other_3",
               n_participants = n_distinct(x$RID),
               pet_positive_last_n = n_distinct(x$RID[x$last_pet == 1L])),
        boot_metrics(x, cols, seed = seed)
      )
    })
  } else {
    map_dfr(ASSAYS, function(from) map_dfr(setdiff(ASSAYS, from), function(to) {
      x <- dat %>% filter(from_assay == from, to_assay == to)
      if (!nrow(x)) return(tibble())
      bind_cols(
        tibble(from_assay = from, to_assay = to,
               n_participants = n_distinct(x$RID),
               pet_positive_last_n = sum(x$last_pet == 1L)),
        boot_metrics(x, cols, seed = seed)
      )
    }))
  }
}
TRANSITIONS <- c("PET_0_to_1", "PET_0_to_0", "PET_1_to_1")
ELIGIBILITIES <- c("all_pet_matched", "first_source_correct")
by_pair <- list()
by_source <- list()
group_counts <- list()
persistent_pair <- list()
persistent_source <- list()
persistent_counts <- list()
for (i in seq_len(nrow(specs))) {
  s <- specs[i, ]
  allx <- pet_filter(pairs, s$window, s$exclude_near_cutoff)
  expected <- coverage$n_both_pet[[i]]
  if (n_distinct(allx$RID) != expected || nrow(allx) != 12L * expected)
    stop("PET eligibility mismatch between participant and paired-assay datasets.")
  for (el in ELIGIBILITIES) {
    eligible <- if (el == "first_source_correct")
      allx %>% filter(source_first_pet_consistent) else allx
    for (tr in TRANSITIONS) {
      x <- eligible %>% filter(pet_transition == tr)
      tag <- paste(s$subset, el, tr, sep = "__")
      group_counts[[tag]] <- x %>% group_by(from_assay) %>% summarise(
        subset = s$subset, baseline_eligibility = el, pet_transition = tr,
        n_participants = n_distinct(RID),
        n_directed_comparisons = n(),
        .groups = "drop"
      )
      if (!nrow(x)) next
      p <- make_table(x, METRICS, seed = BOOT_SEED + i)
      q <- make_table(x, METRICS, by_source = TRUE, seed = BOOT_SEED + i)
      by_pair[[tag]] <- p %>% mutate(subset = s$subset,
        baseline_eligibility = el, pet_transition = tr, .before = 1L)
      by_source[[tag]] <- q %>% mutate(subset = s$subset,
        baseline_eligibility = el, pet_transition = tr, .before = 1L)
    }
  }
  for (tr in c("ALL", TRANSITIONS)) {
    y <- if (tr == "ALL") allx else allx %>% filter(pet_transition == tr)
    if (!nrow(y)) next
    tag <- paste(s$subset, tr, sep = "__")
    z1 <- make_table(y, PERSISTENCE_METRICS, seed = BOOT_SEED + 20L + i)
    z2 <- make_table(y, PERSISTENCE_METRICS, by_source = TRUE,
                     seed = BOOT_SEED + 20L + i)
    persistent_pair[[tag]] <- z1 %>% mutate(subset = s$subset,
        pet_transition = tr, .before = 1L)
    persistent_source[[tag]] <- z2 %>% mutate(subset = s$subset,
        pet_transition = tr, .before = 1L)
    retention_counts <- y %>% group_by(from_assay, to_assay) %>% summarise(
      n_participants = n_distinct(RID),
      informative_first_n = sum(informative_first),
      informative_last_n = sum(informative_last),
      informative_both_n = sum(informative_both),
      same_favored_n = sum(same_favored_when_informative),
      sign_retention_per_100_informative_both = if_else(
        informative_both_n > 0L,
        100 * same_favored_n / informative_both_n, NA_real_),
      .groups = "drop"
    )
    persistent_counts[[tag]] <- retention_counts %>% mutate(
      subset = s$subset, pet_transition = tr, .before = 1L)
  }
}
bind_or_empty <- function(x) if (length(x)) bind_rows(x) else tibble()
write_out(bind_or_empty(group_counts), "34_03_stratum_denominators_by_source.csv")
write_out(bind_or_empty(by_pair), "34_04_PET_conversion_by_directed_pair.csv")
write_out(bind_or_empty(by_source), "34_05_PET_conversion_by_source_assay.csv")
write_out(bind_or_empty(persistent_pair), "34_06_first_last_assay_advantage_pairwise.csv")
write_out(bind_or_empty(persistent_source), "34_07_first_last_assay_advantage_by_source.csv")
write_out(bind_or_empty(persistent_counts), "34_08_persistence_informative_denominators.csv")
headline_names <- c(
  A_last_positive = "Keep_A_last_positive_per_100",
  B_last_positive = "Switch_B_last_positive_per_100",
  positive_rd = "Switch_minus_keep_last_positive_pp",
  A_apparent_up = "Keep_A_apparent_0to1_per_100",
  B_apparent_up = "Switch_B_apparent_0to1_per_100",
  apparent_up_rd = "Switch_minus_keep_apparent_0to1_pp",
  A_apparent_down = "Keep_A_apparent_1to0_per_100",
  B_apparent_down = "Switch_B_apparent_1to0_per_100",
  apparent_down_rd = "Switch_minus_keep_apparent_1to0_pp",
  last_error_rd = "Switch_minus_keep_PET_mismatch_pp",
  false_definitive_rd = "Switch_minus_keep_false_definitive_pp"
)
headlines <- bind_or_empty(by_source) %>%
  filter(subset == "PET90", baseline_eligibility == "first_source_correct",
         metric %in% names(headline_names)) %>%
  mutate(metric = unname(headline_names[metric])) %>%
  select(subset, baseline_eligibility, pet_transition, from_assay,
         n_participants, metric, per_100, ci_low, ci_high)
write_out(headlines, "34_09_PRIMARY_90day_first_correct_headlines.csv")
persistence_headlines <- bind_or_empty(persistent_source) %>%
  filter(subset == "PET90", pet_transition == "ALL",
         metric %in% c("first_error_rd", "last_error_rd", "delta_error_rd")) %>%
  select(subset, pet_transition, from_assay, n_participants,
         metric, per_100, ci_low, ci_high)
write_out(persistence_headlines, "34_10_PRIMARY_first_last_PET_advantage.csv")
readme <- c(
  "34 - Two-timepoint PET-referenced exploratory analysis",
  paste("Core cache:", file.path(results_dir, "02_core_results.rds")),
  paste("Bootstrap replicates:", N_BOOT, "(participant resampling, fixed OOF thresholds)"),
  "PET 0->1: PET-negative first, PET-positive last (not a proof of clinical conversion).",
  "PET 0->0: PET-negative first and last.",
  "PET 1->1: PET-positive first and last.",
  "PET 1->0: report descriptive audit only; likely too few participants.",
  "For first_source_correct, A FIRST plasma classification matches FIRST PET.",
  "A->A: keep the source platform at LAST; A->B: replace it at LAST.",
  "Both trajectories share A FIRST and the identical LAST blood specimen/PET.",
  "PET90 is the prespecified stricter analysis; PET180 is sensitivity.",
  "PET90_no_CL15to25 excludes anyone with either PET Centiloid in [15,25].",
  "For PET 0->1 under first_source_correct, B_last_positive is detection; 100-minus is missed.",
  "For PET 0->0 under first_source_correct, B_apparent_up is spurious apparent blood conversion.",
  "For PET 1->1 under first_source_correct, B_apparent_down is apparent loss of blood positivity.",
  "Positive differences are B minus A, not proof that switching generally harms accuracy.",
  "Source summaries average the three destinations within a participant before bootstrap.",
  "All 12 reciprocal pair net error differences cancel only in an UNRESTRICTED symmetric set.",
  "After conditioning on source-first PET consistency, reciprocity no longer holds.",
  "PET-referenced advantage uses binary Youden predictions; 3-level false-definitive is secondary.",
  "Sign retention requires DIFFERENT PET-referenced error status between A and B at BOTH times.",
  "If informative_both_n is small, sign retention is unstable; inspect counts before interpretation.",
  "These are exploratory, descriptive conditional bootstrap intervals: no cutoff refitting,",
  "no multiplicity correction, and no real-world observed clinical platform switching.",
  "Do not describe the PET binary transition as validation of biomarker CHANGE DIRECTION.",
  ""
)
writeLines(readme, con = file.path(out_dir, "34_README_interpretation.txt"))
cat("\n=============== BOTH-TIMEPOINT PET COVERAGE ===============\n")
print(coverage, n = Inf, width = Inf)
cat("\n========== PRIMARY CLINICAL METRICS (PET90) ==========\n")
print(headlines %>% filter(metric %in% c(
  "Switch_minus_keep_last_positive_pp",
  "Switch_minus_keep_apparent_0to1_pp",
  "Switch_minus_keep_apparent_1to0_pp",
  "Switch_minus_keep_PET_mismatch_pp"
)), n = Inf, width = Inf)
cat("\n========== FIRST vs LAST RELATIVE PET ADVANTAGE ==========\n")
print(persistence_headlines, n = Inf, width = Inf)
cat("\nDone. New results written to: ", out_dir, "\n", sep = "")
cat("IMPORTANT: Simulated platform replacement; exploratory PET-referenced classification only.\n")
