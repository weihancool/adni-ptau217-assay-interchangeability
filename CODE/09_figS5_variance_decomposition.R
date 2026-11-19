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
if (!exists("study_b_code_dir", inherits = TRUE))
  stop("Run 01_config.R manually first, in this same R session.", call. = FALSE)
if (!sb_have("sb_decision_long"))
  stop("Run 02_preprocess.R and 03_core_computation.R manually first.", call. = FALSE)
sb_assert_code_current()
sb_margin_scales <- c("margin_triage", "margin_confirmatory", "margin_youden")
sb_all_scales <- c(sb_margin_scales, "calibrated_probability")
study_b_assay_as_fixed <- TRUE
sb_fit_vc <- function(data, outcome, audit,
                      assay_as_fixed = study_b_assay_as_fixed,
                      include_assay_fold = FALSE, include_assay_run = FALSE,
                      return_model = FALSE) {
  if (include_assay_run && !isTRUE(get0("study_b_has_run_field", ifnotfound = FALSE))) {
    stop("include_assay_run = TRUE but RUN_NO is absent from the input data.",
         call. = FALSE)
  }
  terms <- c("(1|RID)", "(1|occ_id)", "(1|RID:assay)")
  if (!assay_as_fixed) terms <- c(terms, "(1|assay)")
  if (include_assay_fold) terms <- c(terms, "(1|assay_fold)")
  if (include_assay_run)  terms <- c(terms, "(1|assay_run)")
  fixed <- sb_vc_fixed_terms(data, assay_as_fixed)
  f <- as.formula(paste0(outcome, " ~ ", paste(fixed, collapse = " + "), " + ",
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
  out <- tibble(audit = audit, scale = outcome, component = names(vals),
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
                ratio_definition = "person_by_assay / shared_occasion_within_person",
                person_assay_to_occasion_plus_residual_ratio =
                  sb_safe_ratio(vals[["person_by_assay"]],
                                vals[["occasion_within_person"]] + vals[["residual"]]))
  if (return_model) return(list(table = out, model = m))
  out
}
.sb_primary_youden_fit <- sb_fit_vc(sb_decision_long, "margin_youden",
                                    "crossfit_main", assay_as_fixed = TRUE,
                                    return_model = TRUE)
sb_m_youden_primary <- .sb_primary_youden_fit$model
sb_vc_main <- bind_rows(
  map_dfr(sb_all_scales, function(s) {
    if (identical(s, "margin_youden")) return(.sb_primary_youden_fit$table)
    sb_fit_vc(sb_decision_long, s, "crossfit_main", assay_as_fixed = TRUE)
  }),
  map_dfr(sb_all_scales,
          ~ sb_fit_vc(sb_decision_long, .x, "assay_random",
                      assay_as_fixed = FALSE)),
  map_dfr(sb_all_scales,
          ~ sb_fit_vc(sb_decision_long, .x, "crossfit_plus_assay_by_fold",
                      assay_as_fixed = TRUE, include_assay_fold = TRUE))
)
sb_vc_fig2_all <- sb_vc_main
sb_write(sb_vc_main, "Figure_S5_variance_components_main.csv")
sb_vc_model_audit <- sb_vc_main %>%
  group_by(audit, scale) %>%
  summarise(singular = first(singular), converged = first(converged),
            boundary_components = paste(component[is.finite(variance) & variance <= 1e-10], collapse = "; "),
            convergence_message = first(convergence_message), .groups = "drop")
sb_write(sb_vc_model_audit, "Figure_S5_model_fit_audit.csv")
if (any(sb_vc_model_audit$singular))
  warning(sum(sb_vc_model_audit$singular), " of ", nrow(sb_vc_model_audit),
          " models have a singular fit; inspect Figure_S5_model_fit_audit.csv. ",
          "Converged boundary fits are retained and their near-zero components reported.", call. = FALSE)
sb_ratio_ci <- tibble(
  statistic = c("person_assay_to_shared_occasion_ratio",
                "person_assay_to_occasion_plus_residual_ratio"),
  scale = "margin_youden", estimate = NA_real_,
  low = NA_real_, high = NA_real_,
  bootstrap_successful = 0L, bootstrap_attempted = study_b_boot_n_vc)
.sb_primary_youden <- sb_vc_main %>%
  filter(audit == "crossfit_main", scale == "margin_youden") %>%
  slice(1)
if (nrow(.sb_primary_youden) == 1L) {
  sb_ratio_ci$estimate <- c(
    .sb_primary_youden$person_assay_to_occasion_ratio,
    .sb_primary_youden$person_assay_to_occasion_plus_residual_ratio
  )
}
.sb_ratio_model <- function(d) {
  id <- if (".boot_RID" %in% names(d)) ".boot_RID" else "RID"
  d2 <- d %>%
    mutate(RID = factor(.data[[id]]),
           occ_id = interaction(RID, sample_id, drop = TRUE),
           assay = factor(assay, levels = study_b_assays),
           phase = droplevels(factor(phase)),
           baseline_amyloid_group = droplevels(factor(baseline_amyloid_group)))
  fixed <- c("years_from_baseline * baseline_amyloid_group", "assay")
  if (nlevels(d2$phase) >= 2) fixed <- c(fixed, "phase")
  if (nlevels(d2$baseline_amyloid_group) < 2) {
    fixed <- c("years_from_baseline", "assay")
    if (nlevels(d2$phase) >= 2) fixed <- c(fixed, "phase")
  }
  f <- as.formula(paste0(
    "margin_youden ~ ", paste(fixed, collapse = " + "),
    " + (1|RID) + (1|occ_id) + (1|RID:assay)"))
  m <- sb_checked_lmer(f, data = d2, REML = TRUE,
                       control = lmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5)))
  v <- as.data.frame(VarCorr(m))
  values <- c(
    person_assay_to_shared_occasion_ratio =
      sb_safe_ratio(v$vcov[v$grp == "RID:assay"][1],
                    v$vcov[v$grp == "occ_id"][1]),
    person_assay_to_occasion_plus_residual_ratio =
      sb_safe_ratio(v$vcov[v$grp == "RID:assay"][1],
                    v$vcov[v$grp == "occ_id"][1] +
                      v$vcov[v$grp == "Residual"][1])
  )
  list(values = values,
       singular = isSingular(m),
       converged = sb_lmer_converged(m))
}
sb_ratio_draws <- tibble(
  draw = integer(), statistic = character(), value = double(),
  singular = logical(), converged = logical(), success = logical(), failure_reason = character()
)
if (study_b_boot_n_vc > 0) {
  message("Bootstrapping both variance ratios with ", study_b_boot_n_vc,
          " participant resamples and one model fit per resample.")
  .sb_clusters <- split(sb_decision_long,
                        as.character(sb_decision_long$RID), drop = TRUE)
  .sb_cluster_names <- names(.sb_clusters)
  .sb_n_clusters <- length(.sb_clusters)
  .sb_progress_every <- max(1L, floor(study_b_boot_n_vc / 20L))
  set.seed(sb_seed("vc_ratio_shared"))
  sb_ratio_draws <- map_dfr(seq_len(study_b_boot_n_vc), function(b) {
    picked <- sample(.sb_cluster_names, .sb_n_clusters, replace = TRUE)
    boot_data <- bind_rows(map2(picked, seq_along(picked), function(id, j) {
      z <- .sb_clusters[[id]]
      z$.boot_RID <- sprintf("boot_%05d_%05d", b, j)
      z
    }))
    fit_error <- NA_character_
    fit <- tryCatch(suppressWarnings(.sb_ratio_model(boot_data)),
                    error = function(e) { fit_error <<- conditionMessage(e); NULL })
    if (b %% .sb_progress_every == 0L || b == study_b_boot_n_vc) {
      message("Variance-ratio bootstrap: ", b, "/", study_b_boot_n_vc)
    }
    if (is.null(fit)) {
      return(tibble(
        draw = b,
        statistic = c("person_assay_to_shared_occasion_ratio",
                      "person_assay_to_occasion_plus_residual_ratio"),
        value = NA_real_, singular = NA, converged = FALSE, success = FALSE,
        failure_reason = fit_error
      ))
    }
    tibble(
      draw = b,
      statistic = names(fit$values),
      value = unname(fit$values),
      singular = fit$singular,
      converged = fit$converged,
      success = isTRUE(fit$converged) & is.finite(value) & value >= 0,
      failure_reason = case_when(
        !isTRUE(fit$converged) ~ "Model did not converge",
        !is.finite(value) ~ "Non-finite ratio (including zero denominator)",
        value < 0 ~ "Negative variance ratio",
        TRUE ~ NA_character_)
    )
  })
  sb_ratio_ci <- sb_ratio_ci %>%
    select(statistic, scale, estimate, bootstrap_attempted) %>%
    left_join(
      sb_ratio_draws %>%
        filter(success) %>%
        group_by(statistic) %>%
        summarise(
          low = quantile(value, 0.025, na.rm = TRUE, names = FALSE),
          high = quantile(value, 0.975, na.rm = TRUE, names = FALSE),
          bootstrap_successful = n(), .groups = "drop"
        ),
      by = "statistic"
    ) %>%
    mutate(bootstrap_successful = coalesce(bootstrap_successful, 0L),
           bootstrap_attempted = study_b_boot_n_vc)
  if (any(sb_ratio_ci$bootstrap_successful < 0.80 * study_b_boot_n_vc)) {
    warning("Fewer than 80% of variance-ratio bootstrap fits succeeded for at least one statistic.",
            call. = FALSE)
  }
}
sb_ratio_ci$CI_scope <- "Participant model-refitting bootstrap conditional on frozen calibration"
sb_write(sb_ratio_ci, "Figure_S5_variance_ratio_ci.csv")
sb_write(sb_ratio_draws, "Figure_S5_variance_ratio_bootstrap_draws.csv")
sb_ratio_bootstrap_audit <- sb_ratio_draws %>% group_by(statistic) %>%
  summarise(attempted = n(), retained = sum(success),
            retained_singular = sum(success & singular, na.rm = TRUE),
            rejected = sum(!success), .groups = "drop")
sb_write(sb_ratio_bootstrap_audit, "Figure_S5_variance_ratio_bootstrap_audit.csv")
sb_write(sb_ratio_draws %>% filter(!success), "Figure_S5_variance_ratio_bootstrap_failures.csv")
print(sb_ratio_ci, n = Inf, width = Inf)
print(sb_ratio_bootstrap_audit, n = Inf, width = Inf)
sb_person_assay_modes <- local({
  re <- lme4::ranef(sb_m_youden_primary, condVar = TRUE)
  if (!"RID:assay" %in% names(re)) {
    stop("The primary Youden fit does not contain a 'RID:assay' grouping ",
         "factor; panel C cannot be drawn.", call. = FALSE)
  }
  tab <- re[["RID:assay"]]
  pv <- attr(tab, "postVar")
  if (is.null(pv)) {
    stop("Conditional variances are unavailable; refit with condVar = TRUE.",
         call. = FALSE)
  }
  se <- sqrt(as.numeric(pv))
  key <- rownames(tab)
  parts <- strsplit(key, ":", fixed = TRUE)
  p1 <- vapply(parts, function(z) z[1], character(1))
  p2 <- vapply(parts, function(z) if (length(z) >= 2) z[2] else NA_character_,
               character(1))
  assay_is_second <- p2 %in% as.character(study_b_assays)
  assay_vec <- ifelse(assay_is_second, p2, p1)
  rid_vec   <- ifelse(assay_is_second, p1, p2)
  if (!all(assay_vec %in% as.character(study_b_assays))) {
    stop("Could not parse assay labels out of the RID:assay random-effect ",
         "level names.", call. = FALSE)
  }
  tibble(
    RID = rid_vec,
    assay = assay_vec,
    conditional_mode = as.numeric(tab[[1]]),
    conditional_se = se
  ) %>%
    mutate(low = conditional_mode - 1.96 * conditional_se,
           high = conditional_mode + 1.96 * conditional_se,
           excludes_zero = (low > 0) | (high < 0)) %>%
    arrange(conditional_mode) %>%
    mutate(rank = row_number())
})
sb_write(sb_person_assay_modes, "Figure_S5_person_by_assay_conditional_modes.csv")
sb_switch_flip <- local({
  wide <- sb_decision_long %>%
    select(RID_num, sample_id, assay, margin_youden) %>%
    mutate(assay = as.character(assay)) %>%
    tidyr::pivot_wider(names_from = assay, values_from = margin_youden)
  absent <- setdiff(as.character(study_b_assays), names(wide))
  if (length(absent)) {
    stop("Panel D is missing assay column(s): ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  ordered_pairs <- do.call(rbind, lapply(sb_pairs, function(p) {
    rbind(data.frame(reference = p[1], other = p[2], stringsAsFactors = FALSE),
          data.frame(reference = p[2], other = p[1], stringsAsFactors = FALSE))
  }))
  long <- map_dfr(seq_len(nrow(ordered_pairs)), function(k) {
    a <- ordered_pairs$reference[k]
    b <- ordered_pairs$other[k]
    tibble(RID = wide$RID_num, sample_id = wide$sample_id,
           reference_assay = a, other_assay = b,
           reference_margin = wide[[a]], other_margin = wide[[b]])
  }) %>%
    filter(is.finite(reference_margin), is.finite(other_margin)) %>%
    mutate(distance = abs(reference_margin),
           flipped = (reference_margin >= 0) != (other_margin >= 0))
  breaks <- c(seq(0, 2, by = 0.25), Inf)
  long %>%
    mutate(bin = cut(distance, breaks = breaks, include.lowest = TRUE,
                     right = FALSE)) %>%
    group_by(bin) %>%
    summarise(
      comparisons = n(),
      participants = n_distinct(RID),
      distance_median = median(distance),
      flip_rate = mean(flipped),
      .groups = "drop"
    ) %>%
    filter(comparisons >= 30) %>%
    mutate(share_of_comparisons = comparisons / sum(comparisons),
           unit = "ordered assay pairs x collection occasions")
})
sb_write(sb_switch_flip, "Figure_S5_assay_switch_flip_by_margin.csv")
.vc_youden <- sb_vc_main %>%
  filter(audit == "crossfit_main", scale == "margin_youden")
.v_pa  <- .vc_youden$variance[.vc_youden$component == "person_by_assay"][1]
.v_res <- .vc_youden$variance[.vc_youden$component == "residual"][1]
sd_switch_total <- if (is.finite(.v_pa) && is.finite(.v_res))
  sqrt(2 * (.v_pa + .v_res)) else NA_real_
sd_switch_person <- if (is.finite(.v_pa)) sqrt(2 * .v_pa) else NA_real_
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required to assemble Figure 2.", call. = FALSE)
}
FIG_W_MM <- 190
FIG_H_MM <- 152
FONT_FAMILY <- "sans"
COL_TEXT <- "#263238"
COL_AXIS <- "#66727A"
component_order <- c(
  "person", "person_by_assay", "occasion_within_person",
  "assay", "assay_by_fold", "assay_by_run", "residual"
)
component_labels <- c(
  person = "Person",
  person_by_assay = "Person × assay",
  occasion_within_person = "Occasion",
  assay = "Assay",
  assay_by_fold = "Assay × fold",
  assay_by_run = "Assay × run",
  residual = "Residual"
)
component_colours <- c(
  person = "#305F78",
  person_by_assay = "#D49A45",
  occasion_within_person = "#8FB4C6",
  assay = "#A6AEB5",
  assay_by_fold = "#8C78A8",
  assay_by_run = "#6F9E88",
  residual = "#E3E7E9"
)
scale_order <- c(
  "margin_triage", "margin_confirmatory", "margin_youden",
  "calibrated_probability"
)
scale_labels <- c(
  margin_triage = "Triage margin",
  margin_confirmatory = "Confirmatory margin",
  margin_youden = "Youden margin",
  calibrated_probability = "Calibrated probability"
)
assay_display <- function(x) {
  dplyr::recode(as.character(x),
                "AlzPath" = "ALZpath", "ALZpath" = "ALZpath",
                "C2N_ratio" = "C2N", "C2N" = "C2N",
                .default = as.character(x))
}
assay_levels <- assay_display(study_b_assays)
assay_colours <- setNames(c("#305F78", "#D49A45", "#7FA9BC", "#9A7FA8"),
                          assay_levels)
theme_fig2 <- sb_theme_classic(base_size = 8.5, base_family = FONT_FAMILY) +
  theme(
    text = element_text(colour = COL_TEXT),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.8, colour = COL_TEXT),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.line = element_line(colour = COL_AXIS, linewidth = 0.3),
    axis.ticks = element_line(colour = COL_AXIS, linewidth = 0.3),
    axis.ticks.length = grid::unit(1.2, "mm"),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(7, 8, 6, 7),
    plot.tag = element_text(size = 11, face = "bold", colour = COL_TEXT),
    plot.tag.position = "topleft",
    legend.title = element_blank(),
    legend.text = element_text(size = 7.4),
    legend.key.height = grid::unit(3.0, "mm"),
    legend.key.width = grid::unit(4.4, "mm")
  )
dat_A <- sb_vc_main %>%
  filter(audit == "crossfit_main", !is.na(variance)) %>%
  mutate(
    scale = factor(scale, levels = rev(scale_order),
                   labels = rev(unname(scale_labels[scale_order]))),
    component = factor(component, levels = component_order),
    label = if_else(proportion >= 0.09,
                    scales::percent(proportion, accuracy = 0.1), ""),
    label_colour = if_else(component %in% c("person"), "white", COL_TEXT)
  ) %>%
  arrange(scale, component) %>%
  group_by(scale) %>%
  mutate(label_x = cumsum(proportion) - proportion / 2) %>%
  ungroup()
pA <- ggplot(dat_A, aes(proportion, scale, fill = component)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.25,
           position = position_stack(reverse = TRUE)) +
  geom_text(aes(x = label_x, label = label, colour = label_colour),
            size = 2.3, family = FONT_FAMILY) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = scales::label_percent(accuracy = 1),
                     expand = expansion(mult = 0)) +
  scale_fill_manual(values = component_colours, breaks = component_order,
                    labels = component_labels, drop = TRUE) +
  scale_colour_identity() +
  labs(x = "Share of total random variance", y = NULL, tag = "A") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_fig2 +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.justification = "left",
        legend.box.just = "left", legend.margin = margin(2, 0, 0, 0))
ratio_labels <- c(
  person_assay_to_shared_occasion_ratio = "Person × assay / occasion",
  person_assay_to_occasion_plus_residual_ratio =
    "Person × assay / (occasion + residual)"
)
ratio_order <- names(ratio_labels)
dat_B <- sb_ratio_ci %>%
  mutate(
    y = length(ratio_order) + 1L - match(statistic, ratio_order),
    ci_available = is.finite(low) & is.finite(high),
    result_label = if_else(ci_available,
                           sprintf("%.2f  [%.2f, %.2f]", estimate, low, high),
                           sprintf("%.2f", estimate))
  ) %>%
  filter(is.finite(estimate))
draw_B <- sb_ratio_draws %>%
  filter(success, statistic %in% ratio_order, is.finite(value)) %>%
  mutate(y = length(ratio_order) + 1L - match(statistic, ratio_order))
draw_summary_B <- draw_B %>%
  group_by(statistic) %>%
  summarise(q25 = quantile(value, 0.25, na.rm = TRUE, names = FALSE),
            q75 = quantile(value, 0.75, na.rm = TRUE, names = FALSE),
            .groups = "drop")
dat_B <- dat_B %>% left_join(draw_summary_B, by = "statistic")
density_B <- bind_rows(lapply(split(draw_B, draw_B$statistic), function(z) {
  if (nrow(z) < 2L || stats::sd(z$value) == 0) return(tibble())
  den <- stats::density(z$value, n = 512, na.rm = TRUE, cut = 0)
  y0 <- unique(z$y)[1]
  tibble(statistic = unique(z$statistic)[1], x = den$x,
         ymin = y0, ymax = y0 + 0.32 * den$y / max(den$y))
}))
if (nrow(density_B) == 0L) {
  density_B <- tibble(statistic = character(), x = double(),
                      ymin = double(), ymax = double())
}
pB <- ggplot() +
  geom_vline(xintercept = 1, colour = "#9AA3A9", linewidth = 0.45,
             linetype = "22") +
  geom_ribbon(data = density_B,
              aes(x = x, ymin = ymin, ymax = ymax, group = statistic,
                  fill = statistic), alpha = 0.55, colour = NA) +
  geom_segment(data = filter(dat_B, ci_available),
               aes(x = low, xend = high, y = y - 0.07, yend = y - 0.07,
                   colour = statistic), linewidth = 0.65, lineend = "round") +
  geom_segment(data = filter(dat_B, is.finite(q25), is.finite(q75)),
               aes(x = q25, xend = q75, y = y - 0.07, yend = y - 0.07,
                   colour = statistic), linewidth = 1.8, lineend = "round") +
  geom_point(data = dat_B, aes(x = estimate, y = y - 0.07, colour = statistic),
             shape = 21, size = 2.6, stroke = 0.75, fill = "white") +
  geom_text(data = dat_B, aes(x = Inf, y = y + 0.26, label = result_label),
            hjust = 1.05, size = 2.4, family = FONT_FAMILY, colour = COL_TEXT) +
  scale_fill_manual(values = c(
    "person_assay_to_shared_occasion_ratio" = "#D49A45",
    "person_assay_to_occasion_plus_residual_ratio" = "#8FB4C6")) +
  scale_colour_manual(values = c(
    "person_assay_to_shared_occasion_ratio" = "#D49A45",
    "person_assay_to_occasion_plus_residual_ratio" = "#305F78")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  scale_y_continuous(breaks = seq_along(ratio_order),
                     labels = rev(unname(ratio_labels[ratio_order])),
                     limits = c(0.62, length(ratio_order) + 0.48),
                     expand = expansion(mult = 0)) +
  labs(x = "Variance ratio (Youden margin)", y = NULL, tag = "B") +
  guides(fill = "none", colour = "none") +
  theme_fig2 +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
dat_C <- sb_person_assay_modes %>%
  mutate(assay_label = factor(assay_display(assay), levels = assay_levels))
pct_excluding_zero <- 100 * mean(dat_C$excludes_zero)
mode_sd <- stats::sd(dat_C$conditional_mode)
n_modes <- nrow(dat_C)
pC <- ggplot(dat_C, aes(x = rank)) +
  geom_hline(yintercept = 0, colour = COL_TEXT, linewidth = 0.35) +
  geom_linerange(aes(ymin = low, ymax = high, colour = assay_label),
                 linewidth = 0.22, alpha = 0.38) +
  geom_point(aes(y = conditional_mode, colour = assay_label),
             size = 0.30, alpha = 0.95) +
  scale_colour_manual(values = assay_colours, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  scale_y_continuous(expand = expansion(mult = 0.04)) +
  annotate("text", x = 1, y = Inf, hjust = -0.03, vjust = 1.5,
           size = 2.4, family = FONT_FAMILY, colour = COL_TEXT,
           label = sprintf(
             "%d participant × assay effects\nSD = %.2f SD units; %.0f%% with 95%% conditional interval excluding 0",
             n_modes, mode_sd, pct_excluding_zero)) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE,
                               override.aes = list(size = 1.6, alpha = 1,
                                                   linewidth = 0.8))) +
  labs(x = "Participant × assay effects, ranked",
       y = "Conditional mode (SD units)", tag = "C") +
  theme_fig2 +
  theme(legend.position = "bottom", legend.justification = "left",
        legend.margin = margin(2, 0, 0, 0))
p <- (pA | pB) / pC + patchwork::plot_layout(heights = c(1, 1.15)) +
  patchwork::plot_annotation(theme = theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(4, 4, 3, 4)))
out_base <- file.path(study_b_dir, "Figure_S5_variance_decomposition")
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
sb_ggsave(paste0(out_base, ".pdf"), p, width = FIG_W_MM, height = FIG_H_MM,
          units = "mm", device = pdf_device, bg = "white")
sb_ggsave(paste0(out_base, ".tiff"), p, width = FIG_W_MM, height = FIG_H_MM,
          units = "mm", dpi = 600, compression = "lzw", bg = "white")
sb_ggsave(paste0(out_base, ".png"), p, width = FIG_W_MM, height = FIG_H_MM,
          units = "mm", dpi = 600, bg = "white")
cat("\nFigure S5 -- variance decomposition (all specifications)\n")
print(sb_vc_main, n = 120)
cat("\nFigure S5 panel C -- participant-by-assay conditional modes\n")
cat(sprintf("  n = %d; SD = %.3f; %% excluding zero = %.1f\n",
            n_modes, mode_sd, pct_excluding_zero))
cat("\nSupporting CSV only -- observed same-draw flips\n")
print(sb_switch_flip)
cat(sprintf("  Switch-shift SD implied by the fit: person x assay %.3f, total %.3f\n",
            sd_switch_person, sd_switch_total))
writeLines(paste0("Figure S5. Supporting variance models. A: fitted random variance shares ",
                  "with assay modelled as a fixed effect; fixed-effect variation is not in the denominator. ",
                  "B: Youden variance ratios and participant model-refitting bootstrap 95% intervals, ",
                  "conditional on frozen calibration. Failed fits and non-finite ratios are excluded; ",
                  "converged singular fits are retained. C: participant-by-assay conditional modes ",
                  "with model-conditional intervals; these are shrunken estimates, not independent ",
                  "participant significance tests or direct proof of longitudinal persistence. ",
                  "The shared occasion component is not a measure of pure biological change. ",
                  "Same-draw margin-dependent reclassification is displayed in Figure 2B."),
           file.path(study_b_dir, "Figure_S5_caption.txt"))
