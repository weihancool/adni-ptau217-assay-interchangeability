suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr) })
root <- Sys.getenv("STUDY_B_ROOT")
res <- if (exists("study_b_dir")) study_b_dir else file.path(root, "RESULTS")
assays <- c("AlzPath", "C2N_ratio", "Fujirebio", "Janssen")
B <- as.integer(Sys.getenv("STUDY_B_UPGRADE_BOOT", unset = "2000"))
seed <- 20260917L
core <- sb_read_core(file.path(res, "02_core_results.rds"))
fl <- core$sb_rep_dat %>%
  mutate(assay = as.character(assay), class = as.character(three_class)) %>%
  arrange(RID, assay, plasma_date, sample_id) %>% group_by(RID, assay) %>%
  summarise(first_margin = first(margin_youden), last_margin = last(margin_youden),
            first_class = first(class), last_class = last(class),
            first_date = first(plasma_date), last_date = last(plasma_date),
            first_sample = first(sample_id), last_sample = last(sample_id), .groups = "drop")
pairs <- map_dfr(assays, function(from) map_dfr(setdiff(assays, from), function(to) {
  s <- fl %>% filter(assay == from) %>%
    rename_with(~paste0("source_", .x), -c(RID, assay)) %>% select(-assay)
  t <- fl %>% filter(assay == to) %>%
    rename_with(~paste0("target_", .x), -c(RID, assay)) %>% select(-assay)
  inner_join(s, t, by = "RID") %>%
    transmute(RID, from_assay = from, to_assay = to,
              source_first_date, target_first_date, source_last_date, target_last_date,
              source_first_sample, target_first_sample, source_last_sample, target_last_sample,
              delta_same = source_last_margin - source_first_margin,
              delta_switch = target_last_margin - source_first_margin,
              delta_platform = target_last_margin - source_last_margin,
              direction_reversal = delta_same * delta_switch < 0,
              platform_exceeds_same_change = abs(delta_platform) > abs(delta_same))
}))
stopifnot(
  all(pairs$source_first_date == pairs$target_first_date),
  all(pairs$source_last_date == pairs$target_last_date),
  all(pairs$source_first_sample == pairs$target_first_sample),
  all(pairs$source_last_sample == pairs$target_last_sample),
  all(pairs$source_first_sample != pairs$source_last_sample),
  all(pairs$source_last_date > pairs$source_first_date)
)
participant_ci <- function(dat, metric, boot_seed) {
  id <- dat %>% group_by(RID) %>%
    summarise(value = mean(as.numeric(.data[[metric]]), na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(value))
  if (!nrow(id)) return(tibble(metric = metric, n_participants = 0L,
    n_comparisons = nrow(dat), estimate = NA_real_, low = NA_real_, high = NA_real_,
    bootstrap_successful = 0L, bootstrap_attempted = B))
  set.seed(boot_seed)
  draws <- replicate(B, mean(sample(id$value, nrow(id), replace = TRUE)))
  tibble(metric = metric, n_participants = nrow(id), n_comparisons = nrow(dat),
         estimate = mean(id$value),
         low = quantile(draws, .025, na.rm = TRUE, names = FALSE),
         high = quantile(draws, .975, na.rm = TRUE, names = FALSE),
         bootstrap_successful = sum(is.finite(draws)), bootstrap_attempted = B)
}
metrics <- c("direction_reversal", "platform_exceeds_same_change")
cutoffs <- c(0, 0.05, 0.10, 0.20)
threshold_out <- map_dfr(seq_along(cutoffs), function(i) {
  cut <- cutoffs[i]
  dat <- if (cut == 0) pairs else filter(pairs, abs(delta_same) >= cut)
  map_dfr(seq_along(metrics), function(j) participant_ci(dat, metrics[j], seed + 6000L + 10L*i + j)) %>%
    mutate(minimum_abs_same_change_SD = cut,
           restriction = if (cut == 0) "None (frozen primary population)"
                         else paste0("Exclude |delta_same| < ", format(cut, nsmall = 2), " SD"),
           .before = 1)
}) %>% mutate(across(c(estimate, low, high), ~100*.x, .names = "per_100_{.col}"))
source_change <- pairs %>% filter(is.finite(delta_same)) %>% distinct(RID, from_assay, delta_same) %>%
  mutate(abs_same_change = abs(delta_same), same_change_quartile = ntile(abs_same_change, 4L))
quartile_pairs <- pairs %>% filter(is.finite(delta_same)) %>%
  left_join(source_change %>% select(RID, from_assay, same_change_quartile),
            by = c("RID", "from_assay"))
quartile_out <- quartile_pairs %>% group_split(same_change_quartile) %>%
  imap_dfr(function(dat, i) {
    q <- unique(dat$same_change_quartile)
    bounds <- range(abs(dat$delta_same), na.rm = TRUE)
    map_dfr(seq_along(metrics), function(j) participant_ci(dat, metrics[j], seed + 7000L + 10L*q + j)) %>%
      mutate(same_change_quartile = q, min_abs_same_change_SD = bounds[1],
             max_abs_same_change_SD = bounds[2], .before = 1)
  }) %>% mutate(across(c(estimate, low, high), ~100*.x, .names = "per_100_{.col}"))
frozen_path <- file.path(res, "Table_S9_A3_longitudinal_pseudo_switch.csv")
if (file.exists(frozen_path)) {
  frozen <- sb_read_result(frozen_path, check.names = FALSE) %>%
    filter(from_assay == "All", to_assay == "All", metric %in% metrics) %>%
    select(metric, frozen_estimate = estimate, frozen_low = low, frozen_high = high,
           frozen_bootstrap_successful = bootstrap_successful,
           frozen_bootstrap_attempted = bootstrap_attempted)
  chk <- threshold_out %>% filter(minimum_abs_same_change_SD == 0) %>%
    select(metric, new_estimate = estimate) %>% inner_join(frozen, by = "metric")
  stopifnot("A3 sensitivity failed to reproduce frozen primary point estimates" =
              nrow(chk) == length(metrics) &&
              all(abs(chk$new_estimate - chk$frozen_estimate) < 1e-12))
  threshold_out <- threshold_out %>%
    left_join(frozen, by = "metric") %>%
    mutate(low = if_else(minimum_abs_same_change_SD == 0, frozen_low, low),
           high = if_else(minimum_abs_same_change_SD == 0, frozen_high, high),
           bootstrap_successful = if_else(minimum_abs_same_change_SD == 0,
                                          frozen_bootstrap_successful, bootstrap_successful),
           bootstrap_attempted = if_else(minimum_abs_same_change_SD == 0,
                                         frozen_bootstrap_attempted, bootstrap_attempted),
           per_100_low = 100 * low, per_100_high = 100 * high) %>%
    select(-starts_with("frozen_"))
}
sb_save_csv(threshold_out, file.path(res, "Table_S9_A3_small_change_threshold_sensitivity.csv"),
          row.names = FALSE, na = "")
sb_save_csv(quartile_out, file.path(res, "Table_S9_A3_same_change_quartile_sensitivity.csv"),
          row.names = FALSE, na = "")
message("Wrote A3 small-change sensitivity tables.")
