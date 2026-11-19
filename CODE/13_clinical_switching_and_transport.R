options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
})
ROOT <- if (exists("study_b_package_root")) study_b_package_root else Sys.getenv("STUDY_B_ROOT")
RESULTS <- if (exists("study_b_dir")) study_b_dir else file.path(ROOT, "RESULTS")
RAW <- if (exists("study_b_data_dir")) study_b_data_dir else
  Sys.getenv("STUDY_B_DATA_DIR", unset = file.path(ROOT, "01_DATA", "raw_rda"))
N_BOOT <- as.integer(Sys.getenv("STUDY_B_UPGRADE_BOOT", unset = "2000"))
B <- N_BOOT
RUN_LEGACY_B2 <- identical(Sys.getenv("STUDY_B_RUN_LEGACY_B2", unset = "FALSE"), "TRUE")
SEED <- 20260917L
ASSAYS <- c("AlzPath", "C2N_ratio", "Fujirebio", "Janssen")
needed <- file.path(RESULTS, c("02_core_results.rds",
                               "Table_1_participant_audit.csv"))
if (any(!file.exists(needed))) stop("Missing frozen Study B result(s): ",
                                    paste(basename(needed[!file.exists(needed)]), collapse = ", "))
core <- readRDS(needed[1])
scored <- core$sb_scored
anchor <- core$sb_anchor
repdat <- core$sb_rep_dat
persons <- read.csv(needed[2], check.names = FALSE)
stopifnot(nrow(scored) > 0, nrow(anchor) > 0,
          all(ASSAYS %in% unique(as.character(scored$assay))))
write_up <- function(x, name) {
  prefix <- case_when(
    startsWith(name, "A1_") ~ "Table_S7_",
    startsWith(name, "A2_") ~ "Table_S8_",
    startsWith(name, "A3_") ~ "Table_S9_",
    startsWith(name, "B1_") | startsWith(name, "B3_") ~ "Table_S10_",
    startsWith(name, "QC_") ~ "AUDIT_",
    startsWith(name, "HEADLINE_") ~ "",
    startsWith(name, "LEGACY_") ~ "",
    TRUE ~ "AUXILIARY_"
  )
  p <- file.path(RESULTS, paste0(prefix, name, ".csv"))
  write.csv(x, p, row.names = FALSE, na = "")
  message("Wrote ", p)
}
qc_events <- core$sb_events %>%
  count(RID, name = "n_occasions") %>%
  count(n_occasions, name = "n_participants") %>%
  mutate(source = "sb_events (candidate blood draws in frozen core)")
qc_scored_occ <- scored %>%
  mutate(assay = as.character(assay)) %>%
  group_by(RID, sample_id) %>%
  summarise(n_assays = n_distinct(assay[is.finite(margin_youden)]), .groups = "drop")
qc_scored <- qc_scored_occ %>%
  count(RID, name = "n_occasions") %>%
  count(n_occasions, name = "n_participants") %>%
  mutate(source = "sb_scored (occasions entering four-assay analyses)")
qc_occasions <- bind_rows(qc_events, qc_scored) %>%
  group_by(source) %>%
  mutate(total_participants = sum(n_participants),
         total_occasions = sum(n_occasions * n_participants)) %>%
  ungroup() %>%
  relocate(source)
write_up(qc_occasions, "QC_occasion_counts")
message("Occasions: candidate = ", nrow(core$sb_events),
        "; scored = ", nrow(qc_scored_occ),
        "; scored with all four finite assays = ", sum(qc_scored_occ$n_assays == 4))
print(qc_occasions)
if (any(qc_scored_occ$n_assays != length(ASSAYS)))
  warning(sum(qc_scored_occ$n_assays != length(ASSAYS)),
          " scored occasion(s) lack four finite assay margins.")
cluster_mean_ci <- function(data, value, groups = character(), B = N_BOOT,
                            seed = SEED) {
  value <- rlang::ensym(value)
  idmeans <- data %>%
    group_by(across(all_of(c(groups, "RID")))) %>%
    summarise(value = mean(as.numeric(!!value), na.rm = TRUE), .groups = "drop")
  chunks <- if (length(groups)) {
    split(idmeans, interaction(idmeans[groups], drop = TRUE))
  } else {
    list(overall = idmeans)
  }
  set.seed(seed)
  map_dfr(chunks, function(d) {
    x <- d$value[is.finite(d$value)]
    draws <- if (length(x)) replicate(B, mean(sample(x, length(x), replace = TRUE))) else NA_real_
    result <- tibble(n_participants = length(x), estimate = mean(x),
                     low = quantile(draws, .025, na.rm = TRUE, names = FALSE),
                     high = quantile(draws, .975, na.rm = TRUE, names = FALSE),
                     bootstrap_successful = sum(is.finite(draws)),
                     bootstrap_attempted = B)
    if (length(groups)) bind_cols(d[1, groups, drop = FALSE], result) else result
  })
}
cluster_pooled_ci <- function(data, value, stratum, B = N_BOOT, seed = SEED) {
  x <- as.numeric(data[[value]])
  ok <- is.finite(x) & !is.na(data[[stratum]])
  rid <- factor(data$RID[ok])
  st <- droplevels(as.factor(data[[stratum]][ok]))
  S <- tapply(x[ok], list(rid, st), sum);    S[is.na(S)] <- 0
  N <- tapply(x[ok], list(rid, st), length); N[is.na(N)] <- 0
  est <- colSums(S) / colSums(N)
  set.seed(seed)
  n_id <- nrow(S)
  draws <- replicate(B, {
    w <- tabulate(sample.int(n_id, n_id, replace = TRUE), nbins = n_id)
    as.numeric(crossprod(w, S)) / as.numeric(crossprod(w, N))
  })
  draws <- matrix(draws, nrow = ncol(S))
  tibble(stratum = colnames(S),
         comparisons = as.integer(colSums(N)),
         participants = as.integer(colSums(N > 0)),
         estimate = as.numeric(est),
         low = apply(draws, 1, quantile, .025, na.rm = TRUE, names = FALSE),
         high = apply(draws, 1, quantile, .975, na.rm = TRUE, names = FALSE),
         bootstrap_successful = as.integer(rowSums(is.finite(draws))),
         bootstrap_attempted = B)
}
make_switches <- function(dat) {
  w <- dat %>%
    transmute(RID, sample_id, plasma_date, amyloid_positive,
              assay = as.character(assay),
              binary = as.integer(youden_positive),
              class = as.character(three_class), margin = margin_youden) %>%
    pivot_wider(names_from = assay, values_from = c(binary, class, margin))
  map_dfr(ASSAYS, function(from) map_dfr(setdiff(ASSAYS, from), function(to) {
    tibble(RID = w$RID, sample_id = w$sample_id, plasma_date = w$plasma_date,
           amyloid_positive = w$amyloid_positive,
           from_assay = from, to_assay = to,
           from_binary = w[[paste0("binary_", from)]],
           to_binary = w[[paste0("binary_", to)]],
           from_class = w[[paste0("class_", from)]],
           to_class = w[[paste0("class_", to)]],
           from_margin = w[[paste0("margin_", from)]],
           to_margin = w[[paste0("margin_", to)]])
  })) %>%
    mutate(binary_flip = from_binary != to_binary,
           three_class_change = from_class != to_class,
           extreme_switch = (from_class == "negative" & to_class == "positive") |
             (from_class == "positive" & to_class == "negative"),
           from_confirmatory = from_class == "intermediate",
           to_confirmatory = to_class == "intermediate",
           confirmatory_change = from_confirmatory != to_confirmatory,
           switch = paste(from_assay, to_assay, sep = " -> "))
}
summarise_switch_metrics <- function(sw, population, scope, seed = SEED) {
  metrics <- c("binary_flip", "three_class_change", "extreme_switch", "confirmatory_change")
  map_dfr(seq_along(metrics), function(k) {
    m <- metrics[k]
    overall <- cluster_mean_ci(sw, !!rlang::sym(m), B = B, seed = seed + k) %>%
      mutate(from_assay = "All", to_assay = "All")
    pairs <- cluster_mean_ci(sw, !!rlang::sym(m), c("from_assay", "to_assay"),
                             B = B, seed = seed + 100L + k)
    bind_rows(overall, pairs) %>%
      mutate(population = population, scope = scope, metric = m,
             per_100 = 100 * estimate, per_100_low = 100 * low,
             per_100_high = 100 * high, .before = 1)
  })
}
switch_anchor <- make_switches(anchor)
switch_all <- make_switches(scored)
a1 <- bind_rows(
  summarise_switch_metrics(switch_anchor, "Full cohort", "Primary anchor"),
  summarise_switch_metrics(switch_all, "Full cohort", "All collection occasions", SEED + 1000L)
) %>%
  mutate(estimand = if_else(from_assay == "All",
                            "Equal-weight average over 12 directed assay switches (participant-cluster bootstrap)",
                            "Single directed assay switch (participant-cluster bootstrap)"))
write_up(a1, "A1_assay_switching_consequences")
margin_breaks <- c(0, .25, .50, .75, 1, 1.25, 1.5, 1.75, 2, Inf)
switch_all_binned <- switch_all %>%
  mutate(source_margin_bin = cut(abs(from_margin), breaks = margin_breaks,
                                 right = FALSE, include.lowest = TRUE))
a1_margin_metrics <- c("binary_flip", "three_class_change", "confirmatory_change")
a1_margin_long <- map_dfr(seq_along(a1_margin_metrics), function(k) {
  cluster_pooled_ci(switch_all_binned, a1_margin_metrics[k], "source_margin_bin",
                    B = N_BOOT, seed = SEED + 1500L + k) %>%
    mutate(metric = a1_margin_metrics[k], .before = 1)
}) %>%
  rename(source_margin_bin = stratum) %>%
  mutate(source_margin_bin = factor(source_margin_bin,
                                    levels = levels(switch_all_binned$source_margin_bin)),
         per_100 = 100 * estimate, per_100_low = 100 * low,
         per_100_high = 100 * high) %>%
  arrange(metric, source_margin_bin)
write_up(a1_margin_long, "A1_switching_risk_by_source_margin_CI")
a1_margin <- a1_margin_long %>%
  select(metric, source_margin_bin, comparisons, participants, estimate, low, high) %>%
  mutate(metric = factor(metric, levels = a1_margin_metrics)) %>%
  arrange(metric) %>%
  pivot_wider(names_from = metric, values_from = c(estimate, low, high),
              names_glue = "{metric}{ifelse(.value == 'estimate', '', paste0('_', .value))}") %>%
  arrange(source_margin_bin)
write_up(a1_margin, "A1_switching_risk_by_source_margin")
vc_file <- file.path(RESULTS, "Figure_S5_variance_components_main.csv")
if (!file.exists(vc_file))
  stop("Missing variance-component source for residual benchmark: ", vc_file,
       "\nRun 09_figS5_variance_decomposition.R first.", call. = FALSE)
vc_main_for_benchmark <- read.csv(vc_file, check.names = FALSE, stringsAsFactors = FALSE)
vc_resid <- vc_main_for_benchmark %>%
  filter(audit == "crossfit_main",
         scale == "margin_youden",
         component == "residual")
if (nrow(vc_resid) != 1L || !is.finite(vc_resid$variance[1]))
  stop("Could not identify exactly one finite primary Youden residual variance.",
       call. = FALSE)
benchmark_residual_variance <- as.numeric(vc_resid$variance[1])
benchmark_pairwise_residual_sd <- sqrt(2 * benchmark_residual_variance)
BENCHMARK_SEED <- 20261007L
benchmark_rows <- switch_all %>%
  filter(is.finite(from_margin), !is.na(binary_flip)) %>%
  mutate(
    source_margin_bin = cut(abs(from_margin), breaks = margin_breaks,
                            right = FALSE, include.lowest = TRUE),
    p_flip_residual_null = pnorm(-abs(from_margin) / benchmark_pairwise_residual_sd),
    excess_row = as.numeric(binary_flip) - p_flip_residual_null
  )
benchmark_point <- benchmark_rows %>%
  group_by(source_margin_bin) %>%
  summarise(
    comparisons = n(),
    participants = n_distinct(RID),
    observed_per100 = 100 * mean(as.numeric(binary_flip)),
    residual_null_per100 = 100 * mean(p_flip_residual_null),
    excess_pp = observed_per100 - residual_null_per100,
    .groups = "drop"
  )
cluster_excess_ci <- function(dat, B = N_BOOT, seed = BENCHMARK_SEED) {
  z <- dat %>%
    group_by(RID, source_margin_bin) %>%
    summarise(S = sum(excess_row), N = n(), .groups = "drop")
  bins <- levels(dat$source_margin_bin)
  ids <- sort(unique(z$RID))
  S <- matrix(0, nrow = length(ids), ncol = length(bins),
              dimnames = list(as.character(ids), bins))
  N <- matrix(0, nrow = length(ids), ncol = length(bins),
              dimnames = list(as.character(ids), bins))
  ii <- match(z$RID, ids)
  jj <- match(as.character(z$source_margin_bin), bins)
  for (k in seq_len(nrow(z))) {
    S[ii[k], jj[k]] <- z$S[k]
    N[ii[k], jj[k]] <- z$N[k]
  }
  point <- colSums(S) / colSums(N)
  set.seed(seed)
  draws <- replicate(B, {
    sampled <- sample.int(length(ids), length(ids), replace = TRUE)
    w <- tabulate(sampled, nbins = length(ids))
    num <- as.numeric(crossprod(w, S))
    den <- as.numeric(crossprod(w, N))
    num / den
  })
  draws <- matrix(draws, nrow = length(bins))
  tibble(
    source_margin_bin = factor(bins, levels = bins),
    excess_pp_boot = 100 * point,
    excess_low = 100 * apply(draws, 1, quantile, .025,
                             na.rm = TRUE, names = FALSE),
    excess_high = 100 * apply(draws, 1, quantile, .975,
                              na.rm = TRUE, names = FALSE),
    bootstrap_successful = as.integer(rowSums(is.finite(draws))),
    bootstrap_attempted = B
  )
}
benchmark_ci <- cluster_excess_ci(benchmark_rows)
benchmark_observed_ref <- a1_margin_long %>%
  filter(metric == "binary_flip") %>%
  transmute(
    source_margin_bin = as.character(source_margin_bin),
    observed_frozen_per100 = per_100,
    observed_low = per_100_low,
    observed_high = per_100_high
  )
benchmark_final <- benchmark_point %>%
  mutate(source_margin_bin = as.character(source_margin_bin)) %>%
  left_join(benchmark_observed_ref, by = "source_margin_bin") %>%
  left_join(
    benchmark_ci %>% mutate(source_margin_bin = as.character(source_margin_bin)),
    by = "source_margin_bin"
  ) %>%
  mutate(
    observed_difference_from_frozen = observed_per100 - observed_frozen_per100,
    residual_variance = benchmark_residual_variance,
    pairwise_residual_sd = benchmark_pairwise_residual_sd,
    benchmark_definition = paste0(
      "Phi(-abs(source margin)/sqrt(2*residual variance)); residual variance = ",
      signif(benchmark_residual_variance, 6)
    )
  ) %>%
  select(
    source_margin_bin, comparisons, participants,
    observed_per100, observed_low, observed_high,
    residual_null_per100, excess_pp, excess_low, excess_high,
    bootstrap_successful, bootstrap_attempted,
    residual_variance, pairwise_residual_sd,
    observed_difference_from_frozen, benchmark_definition
  )
if (any(!is.finite(benchmark_final$observed_difference_from_frozen)) ||
    max(abs(benchmark_final$observed_difference_from_frozen)) > 1e-10)
  stop("Residual-benchmark observed rates do not reproduce the frozen source-margin estimates.",
       call. = FALSE)
if (any(benchmark_final$bootstrap_successful < .90 * benchmark_final$bootstrap_attempted))
  warning("Residual-benchmark excess bootstrap success below 90% in at least one margin bin.",
          call. = FALSE)
write_up(benchmark_final, "A1_residual_variance_benchmark")
benchmark_output_file <- file.path(RESULTS, "Table_S7_A1_residual_variance_benchmark.csv")
if (exists("sb_stamp_result", mode = "function"))
  sb_stamp_result(benchmark_output_file, extra_inputs = vc_file)
message("Residual-variance benchmark written for Table S7 Panel D / Figure 2B.")
symptomatic_ids <- persons %>%
  filter(diagnosis %in% c("MCI", "Dementia")) %>% pull(RID) %>% as.integer()
if (!length(symptomatic_ids)) stop("No MCI/dementia participants found in audited Table 1 data.")
sym_anchor <- anchor %>% filter(RID %in% symptomatic_ids)
sym_switch <- make_switches(sym_anchor)
a2_switch <- summarise_switch_metrics(sym_switch,
                                      "MCI + dementia (frozen full-cohort OOF thresholds)",
                                      "Primary anchor", SEED + 2000L)
write_up(a2_switch, "A2_symptomatic_switching_consequences")
a2_disagreement <- sym_anchor %>%
  group_by(RID) %>%
  summarise(binary_disagreement = n_distinct(youden_positive) > 1,
            three_class_disagreement = n_distinct(as.character(three_class)) > 1,
            extreme_conflict = any(as.character(three_class) == "negative") &
              any(as.character(three_class) == "positive"),
            .groups = "drop")
a2_disagreement <- map_dfr(c("binary_disagreement", "three_class_disagreement", "extreme_conflict"),
                           function(m) cluster_mean_ci(a2_disagreement, !!rlang::sym(m),
                                                       B = B, seed = SEED + match(m, names(a2_disagreement))) %>%
                             mutate(metric = m, .before = 1))
write_up(a2_disagreement, "A2_symptomatic_disagreement")
a2_pet <- sym_anchor %>% filter(!is.na(amyloid_positive)) %>%
  mutate(assay = as.character(assay)) %>% group_by(assay) %>%
  summarise(n = n(), amyloid_positive_n = sum(amyloid_positive == 1),
            sensitivity = mean(youden_positive[amyloid_positive == 1] == 1, na.rm = TRUE),
            specificity = mean(youden_positive[amyloid_positive == 0] == 0, na.rm = TRUE),
            false_negative_rate = mean(youden_positive[amyloid_positive == 1] == 0, na.rm = TRUE),
            false_positive_rate = mean(youden_positive[amyloid_positive == 0] == 1, na.rm = TRUE),
            .groups = "drop")
write_up(a2_pet, "A2_symptomatic_PET_performance")
a2_vc <- tibble(note = "lme4 unavailable or model not estimable")
if (requireNamespace("lme4", quietly = TRUE)) {
  d <- core$sb_decision_long %>% filter(RID %in% symptomatic_ids) %>%
    mutate(RID = factor(RID), assay = factor(assay, levels = ASSAYS),
           phase = factor(phase), occ_id = interaction(RID, sample_id, drop = TRUE),
           baseline_amyloid_group = factor(baseline_amyloid_group))
  fit <- lme4::lmer(margin_youden ~ assay + years_from_baseline * baseline_amyloid_group +
                      phase + (1|RID) + (1|occ_id) + (1|RID:assay),
                    data = d, REML = TRUE,
                    control = lme4::lmerControl(optimizer = "bobyqa",
                                                optCtrl = list(maxfun = 2e5)))
  vc <- as.data.frame(lme4::VarCorr(fit))
  a2_vc <- vc %>% transmute(component = grp, variance = vcov,
                            proportion = vcov / sum(vcov),
                            n_participants = n_distinct(d$RID),
                            n_occasions = n_distinct(d$occ_id),
                            singular = lme4::isSingular(fit),
                            convergence_message = paste(fit@optinfo$conv$lme4$messages,
                                                        collapse = "; "))
}
write_up(a2_vc, "A2_symptomatic_variance_components")
fl <- repdat %>% mutate(assay = as.character(assay),
                        class = as.character(three_class),
                        binary = as.integer(youden_positive)) %>%
  arrange(RID, assay, plasma_date, sample_id) %>% group_by(RID, assay) %>%
  summarise(first_margin = first(margin_youden), last_margin = last(margin_youden),
            first_class = first(class), last_class = last(class),
            first_binary = first(binary), last_binary = last(binary),
            first_date = first(plasma_date), last_date = last(plasma_date),
            first_sample = first(sample_id), last_sample = last(sample_id), .groups = "drop")
stopifnot(!anyNA(fl$first_date), !anyNA(fl$last_date))
a3_pairs <- map_dfr(ASSAYS, function(from) map_dfr(setdiff(ASSAYS, from), function(to) {
  s <- fl %>% filter(assay == from) %>% rename_with(~paste0("source_", .x),
                                                    -c(RID, assay)) %>% select(-assay)
  t <- fl %>% filter(assay == to) %>% rename_with(~paste0("target_", .x),
                                                  -c(RID, assay)) %>% select(-assay)
  inner_join(s, t, by = "RID") %>%
    transmute(RID, from_assay = from, to_assay = to,
              source_first_date, target_first_date,
              source_last_date, target_last_date,
              source_first_sample, target_first_sample,
              source_last_sample, target_last_sample,
              delta_same = source_last_margin - source_first_margin,
              delta_switch = target_last_margin - source_first_margin,
              delta_platform = target_last_margin - source_last_margin,
              direction_reversal = delta_same * delta_switch < 0,
              direction_reversal_legacy_signdiff = sign(delta_switch) != sign(delta_same),
              platform_exceeds_same_change = abs(delta_platform) > abs(delta_same),
              same_transition = paste(source_first_class, source_last_class, sep = " -> "),
              switched_transition = paste(source_first_class, target_last_class, sep = " -> "),
              same_assay_three_class_changed = source_first_class != source_last_class,
              switched_three_class_changed = source_first_class != target_last_class,
              switch_alters_three_class_transition = same_transition != switched_transition,
              category_change_excess_due_to_switch =
                as.numeric(switched_three_class_changed) -
                as.numeric(same_assay_three_class_changed),
              followup_years = as.numeric(source_last_date - source_first_date) / 365.25)
}))
a3_date_audit <- a3_pairs %>%
  group_by(from_assay, to_assay) %>%
  summarise(n = n(),
            first_date_match = sum(source_first_date == target_first_date),
            last_date_match = sum(source_last_date == target_last_date),
            pct_first_date_match = 100 * mean(source_first_date == target_first_date),
            pct_last_date_match = 100 * mean(source_last_date == target_last_date),
            pct_first_sample_match = 100 * mean(source_first_sample == target_first_sample),
            pct_last_sample_match = 100 * mean(source_last_sample == target_last_sample),
            n_single_draw = sum(source_first_sample == source_last_sample),
            .groups = "drop")
write_up(a3_date_audit, "A3_date_alignment_audit")
stopifnot(
  "A3: 12 directed pairs expected" = nrow(a3_date_audit) == 12L,
  "A3: source/target FIRST blood-draw dates differ" =
    all(a3_pairs$source_first_date == a3_pairs$target_first_date),
  "A3: source/target LAST blood-draw dates differ" =
    all(a3_pairs$source_last_date == a3_pairs$target_last_date),
  "A3: source/target FIRST sample_id differ" =
    all(a3_pairs$source_first_sample == a3_pairs$target_first_sample),
  "A3: source/target LAST sample_id differ" =
    all(a3_pairs$source_last_sample == a3_pairs$target_last_sample)
)
if (any(a3_pairs$source_last_date <= a3_pairs$source_first_date))
  warning("A3: ", n_distinct(a3_pairs$RID[a3_pairs$source_last_date <= a3_pairs$source_first_date]),
          " participant(s) have first date == last date (no follow-up interval).")
a3_zero_audit <- a3_pairs %>%
  group_by(from_assay, to_assay) %>%
  summarise(n = n(),
            n_delta_same_zero = sum(delta_same == 0, na.rm = TRUE),
            n_delta_switch_zero = sum(delta_switch == 0, na.rm = TRUE),
            n_delta_na = sum(is.na(delta_same) | is.na(delta_switch)),
            direction_reversal_strict = mean(direction_reversal, na.rm = TRUE),
            direction_reversal_legacy = mean(direction_reversal_legacy_signdiff, na.rm = TRUE),
            n_definition_discordant = sum(direction_reversal != direction_reversal_legacy_signdiff,
                                          na.rm = TRUE),
            .groups = "drop")
a3_zero_audit <- bind_rows(
  a3_zero_audit,
  a3_zero_audit %>% summarise(from_assay = "All", to_assay = "All",
                              across(c(n, starts_with("n_")), sum),
                              direction_reversal_strict = NA_real_,
                              direction_reversal_legacy = NA_real_))
write_up(a3_zero_audit, "A3_direction_reversal_zero_change_audit")
message("A3 zero-change audit: delta_same == 0: ", sum(a3_pairs$delta_same == 0, na.rm = TRUE),
        "; delta_switch == 0: ", sum(a3_pairs$delta_switch == 0, na.rm = TRUE),
        "; rows where strict and legacy definitions disagree: ",
        sum(a3_pairs$direction_reversal != a3_pairs$direction_reversal_legacy_signdiff, na.rm = TRUE))
a3 <- a3_pairs %>%
  select(-source_first_date, -target_first_date, -source_last_date, -target_last_date,
         -source_first_sample, -target_first_sample, -source_last_sample, -target_last_sample)
a3_metric_names <- c(
  "direction_reversal",
  "platform_exceeds_same_change",
  "same_assay_three_class_changed",
  "switched_three_class_changed",
  "category_change_excess_due_to_switch",
  "switch_alters_three_class_transition"
)
a3_metrics <- map_dfr(a3_metric_names, function(m) {
  overall <- cluster_mean_ci(a3, !!rlang::sym(m), B = B,
                             seed = SEED + 3000L + match(m, names(a3))) %>%
    mutate(from_assay = "All", to_assay = "All")
  pairs <- cluster_mean_ci(a3, !!rlang::sym(m), c("from_assay", "to_assay"), B = B,
                           seed = SEED + 3100L + match(m, names(a3)))
  bind_rows(overall, pairs) %>% mutate(metric = m, per_100 = 100 * estimate,
                                       per_100_low = 100 * low,
                                       per_100_high = 100 * high, .before = 1)
})
write_up(a3_metrics, "A3_longitudinal_pseudo_switch")
write_up(a3 %>% group_by(from_assay, to_assay) %>%
           summarise(n = n(), median_followup_years = median(followup_years),
                     median_abs_same_change = median(abs(delta_same)),
                     median_abs_platform_shift = median(abs(delta_platform)), .groups = "drop"),
         "A3_longitudinal_pseudo_switch_magnitudes")
b1 <- switch_anchor %>% filter(!is.na(amyloid_positive)) %>%
  mutate(binary_direction = case_when(
    !binary_flip ~ "no binary change",
    to_binary == amyloid_positive & from_binary != amyloid_positive ~ "toward PET reference",
    from_binary == amyloid_positive & to_binary != amyloid_positive ~ "away from PET reference",
    TRUE ~ "other"),
    from_rank = match(from_class, c("negative", "intermediate", "positive")),
    to_rank = match(to_class, c("negative", "intermediate", "positive")),
    desired_direction = if_else(amyloid_positive == 1L, 1, -1),
    ordinal_direction = case_when(
      to_rank == from_rank ~ "no three-class change",
      sign(to_rank - from_rank) == desired_direction ~ "toward PET reference",
      TRUE ~ "away from PET reference"))
b1_binary <- b1 %>% count(from_assay, to_assay, amyloid_positive, binary_direction) %>%
  group_by(from_assay, to_assay, amyloid_positive) %>% mutate(proportion = n / sum(n)) %>% ungroup()
b1_ordinal <- b1 %>% count(from_assay, to_assay, amyloid_positive, ordinal_direction) %>%
  group_by(from_assay, to_assay, amyloid_positive) %>% mutate(proportion = n / sum(n)) %>% ungroup()
write_up(b1_binary, "B1_PET_direction_binary_switch")
write_up(b1_ordinal, "B1_PET_direction_three_class_switch")
prevs <- c(.20, .40, .60)
cond_assay <- anchor %>% filter(!is.na(amyloid_positive)) %>%
  mutate(assay = as.character(assay), intermediate = as.character(three_class) == "intermediate") %>%
  group_by(assay, amyloid_positive) %>%
  summarise(p_positive = mean(youden_positive == 1),
            p_negative = mean(as.character(three_class) == "negative"),
            p_intermediate = mean(intermediate),
            p_definitive_positive = mean(as.character(three_class) == "positive"), .groups = "drop")
cond_assay_pos <- cond_assay %>% filter(amyloid_positive == 1) %>%
  select(-amyloid_positive) %>% rename_with(~paste0(.x, "_petpos"), starts_with("p_"))
cond_assay_neg <- cond_assay %>% filter(amyloid_positive == 0) %>%
  select(-amyloid_positive) %>% rename_with(~paste0(.x, "_petneg"), starts_with("p_"))
b3_assay <- crossing(prevalence = prevs, assay = ASSAYS) %>%
  left_join(cond_assay_pos, by = "assay") %>%
  left_join(cond_assay_neg, by = "assay") %>%
  mutate(per100_binary_positive = 100 * (prevalence * p_positive_petpos +
                                           (1-prevalence) * p_positive_petneg),
         per100_false_negative = 100 * prevalence * (1-p_positive_petpos),
         per100_false_positive = 100 * (1-prevalence) * p_positive_petneg,
         per100_definitive_negative = 100 * (prevalence*p_negative_petpos +
                                               (1-prevalence)*p_negative_petneg),
         per100_intermediate_confirmatory = 100 * (prevalence*p_intermediate_petpos +
                                                     (1-prevalence)*p_intermediate_petneg),
         per100_definitive_positive = 100 * (prevalence*p_definitive_positive_petpos +
                                               (1-prevalence)*p_definitive_positive_petneg))
write_up(b3_assay, "B3_hypothetical_prevalence_by_assay")
cond_switch <- b1 %>% group_by(from_assay, to_assay, amyloid_positive) %>%
  summarise(p_binary_flip = mean(binary_flip), p_confirmatory_change = mean(confirmatory_change),
            p_toward_pet = mean(binary_direction == "toward PET reference"),
            p_away_pet = mean(binary_direction == "away from PET reference"), .groups = "drop")
cond_switch_pos <- cond_switch %>% filter(amyloid_positive == 1) %>%
  select(-amyloid_positive) %>% rename_with(~paste0(.x, "_petpos"), starts_with("p_"))
cond_switch_neg <- cond_switch %>% filter(amyloid_positive == 0) %>%
  select(-amyloid_positive) %>% rename_with(~paste0(.x, "_petneg"), starts_with("p_"))
b3_switch <- crossing(prevalence = prevs, from_assay = ASSAYS, to_assay = ASSAYS) %>%
  filter(from_assay != to_assay) %>%
  left_join(cond_switch_pos, by = c("from_assay", "to_assay")) %>%
  left_join(cond_switch_neg, by = c("from_assay", "to_assay")) %>%
  mutate(per100_binary_flip = 100*(prevalence*p_binary_flip_petpos +
                                     (1-prevalence)*p_binary_flip_petneg),
         per100_confirmatory_change = 100*(prevalence*p_confirmatory_change_petpos +
                                             (1-prevalence)*p_confirmatory_change_petneg),
         per100_toward_pet = 100*(prevalence*p_toward_pet_petpos +
                                    (1-prevalence)*p_toward_pet_petneg),
         per100_away_pet = 100*(prevalence*p_away_pet_petpos +
                                  (1-prevalence)*p_away_pet_petneg))
write_up(b3_switch, "B3_hypothetical_prevalence_switching")
stale_b2 <- file.path(RESULTS, c("LEGACY_B2_modifier_data_feasibility.csv",
                                 "LEGACY_B2_BMI_variance_explanation.csv"))
if (any(file.exists(stale_b2)))
  warning("Stale legacy B2 output(s) from an earlier run found; do not report: ",
          paste(basename(stale_b2[file.exists(stale_b2)]), collapse = ", "))
if (RUN_LEGACY_B2) {
  load(file.path(RAW, "VITALS.rda"))
  v <- VITALS %>% transmute(RID = as.integer(RID), date = as.Date(VISDATE),
                            weight = suppressWarnings(as.numeric(VSWEIGHT)),
                            weight_unit = tolower(trimws(as.character(VSWTUNIT))),
                            height = suppressWarnings(as.numeric(VSHEIGHT)),
                            height_unit = tolower(trimws(as.character(VSHTUNIT)))) %>%
    mutate(weight_kg = case_when(grepl("pound", weight_unit) ~ weight*0.45359237,
                                 grepl("kilo", weight_unit) ~ weight, TRUE ~ NA_real_),
           height_m = case_when(grepl("inch", height_unit) ~ height*0.0254,
                                grepl("cent", height_unit) ~ height/100, TRUE ~ NA_real_),
           bmi = weight_kg / height_m^2) %>% filter(is.finite(bmi), bmi >= 10, bmi <= 80)
  base_dates <- persons %>% transmute(RID = as.integer(RID), plasma_date = as.Date(plasma_date))
  bmi_match <- inner_join(base_dates, v, by = "RID", relationship = "many-to-many") %>%
    mutate(gap = abs(as.numeric(date-plasma_date))) %>% arrange(RID, gap) %>% distinct(RID, .keep_all = TRUE)
  load(file.path(RAW, "LABDATA.rda"))
  b2 <- tibble(
    modifier = c("BMI", "eGFR/creatinine", "hemoglobin/anemia"),
    confidently_identified = c(TRUE, TRUE, TRUE),
    participants_within_365_days = c(sum(bmi_match$gap <= 365), NA, NA),
    participants_total_nearest = c(nrow(bmi_match), NA, NA),
    total_anchor_participants = nrow(persons),
    recommendation = c("Feasible as a secondary modifier analysis",
                       "Feasible; run the detailed RCT392/eGFR audit",
                       "Feasible; run the detailed HMT40 audit"),
    reason = c("VITALS contains explicit weight/height fields and units",
               "ADNI/ADNIMERGE mapping: RCT392 = creatinine (rate blanked)",
               "ADNI/ADNIMERGE mapping: HMT40 = hemoglobin"))
  write_up(b2, "LEGACY_B2_modifier_data_feasibility")
  b2_bmi <- tibble(note = "BMI model not fitted")
  if (requireNamespace("lme4", quietly = TRUE)) {
    bmi_cc <- bmi_match %>% filter(gap <= 365) %>%
      transmute(RID, bmi, bmi_z = as.numeric(scale(bmi)))
    db <- core$sb_decision_long %>%
      mutate(RID = as.integer(as.character(RID))) %>%
      inner_join(bmi_cc, by = "RID") %>%
      mutate(RID = factor(RID), assay = factor(assay, levels = ASSAYS),
             phase = factor(phase), occ_id = interaction(RID, sample_id, drop = TRUE),
             baseline_amyloid_group = factor(baseline_amyloid_group))
    fit_bmi <- function(add_interaction) {
      rhs <- paste(
        "assay + years_from_baseline * baseline_amyloid_group + phase + bmi_z",
        if (add_interaction) "+ assay:bmi_z" else "",
        "+ (1|RID) + (1|occ_id) + (1|RID:assay)")
      m <- lme4::lmer(as.formula(paste("margin_youden ~", rhs)), data = db,
                      REML = TRUE,
                      control = lme4::lmerControl(optimizer = "bobyqa",
                                                  optCtrl = list(maxfun = 2e5)))
      vc <- as.data.frame(lme4::VarCorr(m))
      getv <- function(g) vc$vcov[match(g, vc$grp)]
      tibble(model = if (add_interaction) "plus_assay_by_BMI" else "same_subset_unadjusted_for_interaction",
             n_participants = n_distinct(db$RID), n_occasions = n_distinct(db$occ_id),
             person_by_assay_variance = getv("RID:assay"),
             shared_occasion_variance = getv("occ_id"), residual_variance = getv("Residual"),
             person_assay_to_occasion_ratio = person_by_assay_variance/shared_occasion_variance,
             singular = lme4::isSingular(m),
             optimizer_code = if (is.null(m@optinfo$conv$opt)) 0L else m@optinfo$conv$opt[1],
             convergence_message = paste(m@optinfo$conv$lme4$messages, collapse = "; "))
    }
    b2_bmi <- bind_rows(fit_bmi(FALSE), fit_bmi(TRUE)) %>%
      mutate(reference_person_by_assay_variance = first(person_by_assay_variance),
             fraction_remaining = person_by_assay_variance/reference_person_by_assay_variance,
             percent_explained = 100*(1-fraction_remaining))
  }
  write_up(b2_bmi, "LEGACY_B2_BMI365_variance_explanation")
}
headline <- bind_rows(
  a1 %>% filter(scope == "Primary anchor", from_assay == "All") %>%
    transmute(section = "A1 full cohort anchor", metric, estimate, low, high),
  a2_switch %>% filter(from_assay == "All") %>%
    transmute(section = "A2 MCI+dementia anchor", metric, estimate, low, high),
  a2_disagreement %>% transmute(section = "A2 MCI+dementia participant", metric,
                                estimate, low, high),
  a3_metrics %>% filter(from_assay == "All") %>%
    transmute(section = "A3 longitudinal pseudo-switch", metric, estimate, low, high)
) %>% mutate(per_100 = 100*estimate, per_100_low = 100*low,
             per_100_high = 100*high,
             note = case_when(
               grepl("^A1|^A2 MCI\\+dementia anchor", section) ~
                 "Equal-weight average over 12 directed assay switches",
               grepl("^A3", section) & metric == "direction_reversal" ~
                 "Strict definition: delta_same * delta_switch < 0",
               TRUE ~ ""))
write_up(headline, "HEADLINE_analysis_upgrades")
cat("\nAnalysis upgrades complete. Original CODE files were not sourced or modified.\n")
cat("Bootstrap replicates:", N_BOOT, "\n")
print(headline)
