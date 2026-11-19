if (!exists("sb_theme_classic",mode = "function"))
  stop("Run the updated 01_config.R first in this same R session.",call. = FALSE)
.sb_find <- function(file) {
  if (exists("study_b_code_dir", inherits = TRUE)) {
    return(file.path(study_b_code_dir, file))
  }
  env <- Sys.getenv("STUDY_B_ROOT", unset = "")
  if (nzchar(env) && file.exists(file.path(env, "CODE", file))) {
    return(file.path(env, "CODE", file))
  }
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in seq_len(6)) {
    p <- file.path(d, "CODE", file)
    if (file.exists(p)) return(p)
    if (identical(dirname(d), d)) break
    d <- dirname(d)
  }
  stop("Cannot locate CODE/", file, call. = FALSE)
}
if (!exists("study_b_code_dir", inherits = TRUE)) {
  stop("Run 01_config.R manually first, in this same R session.", call. = FALSE)
}
if (!identical(get0("sb_core_contract_used", ifnotfound = NULL), sb_make_core_contract()))
  stop("Figure S6 requires the current 03_core_computation.R objects.", call. = FALSE)
.fig3_reuse <- identical(Sys.getenv("STUDY_B_FIG3_REUSE_OBJECTS", unset = ""), "1") &&
  exists("sb_rep_dat", inherits = TRUE)
if (.fig3_reuse) {
  message("Figure S6: REUSING the sb_rep_dat already in this session ",
          "(STUDY_B_FIG3_REUSE_OBJECTS=1). 03_core_computation.R is NOT ",
          "rerun. Only do this in a session you know is clean.")
} else {
  stop("Run 03_core_computation.R manually first in this same R session; ",
       "automatic sourcing is disabled.", call. = FALSE)
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required to assemble Figure S6.", call. = FALSE)
}
if (!"assay_deviation_adjusted" %in% names(sb_rep_dat)) {
  stop(
    "sb_rep_dat has no assay_deviation_adjusted column.\n",
    "Figure S6 will NOT refit the deviation model locally: a local refit ",
    "silently yields different deviations from the ones behind Table S9 supporting analysis and ",
    "Table S7, which is the discrepancy this revision exists to prevent.\n",
    "Re-run 03_core_computation.R from a clean R session so that it emits ",
    "assay_deviation_adjusted itself.",
    call. = FALSE
  )
}
required_rep_columns <- c(
  "RID", "assay", "plasma_date", "sample_id", "three_class",
  "assay_deviation_adjusted"
)
missing_rep_columns <- setdiff(required_rep_columns, names(sb_rep_dat))
if (length(missing_rep_columns)) {
  stop("sb_rep_dat is missing: ",
       paste(missing_rep_columns, collapse = ", "), call. = FALSE)
}
.dup_occasions <- sb_rep_dat %>%
  count(RID, assay, plasma_date, name = "n_rows") %>%
  filter(n_rows > 1L)
if (nrow(.dup_occasions)) {
  message("Figure S6: ", nrow(.dup_occasions), " RID x assay x plasma_date ",
          "combinations have more than one row; first/last ordering is broken ",
          "by sample_id. Confirm this is expected (for example two protocols ",
          "on the same date).")
} else {
  message("Figure S6: no duplicate RID x assay x plasma_date rows; ",
          "first/last ordering is unambiguous.")
}
sb_dev_wide <- sb_rep_dat %>%
  arrange(RID, assay, plasma_date, sample_id) %>%
  group_by(RID, assay) %>%
  filter(n() >= 2) %>%
  summarise(
    first_deviation = first(assay_deviation_adjusted),
    last_deviation = last(assay_deviation_adjusted),
    first_class = first(as.character(three_class)),
    last_class = last(as.character(three_class)),
    .groups = "drop"
  )
if (nrow(sb_dev_wide) == 0L) {
  stop("No RID x assay cell has at least two observations.", call. = FALSE)
}
.n_zero_dev <- sum(sb_dev_wide$first_deviation == 0, na.rm = TRUE) +
  sum(sb_dev_wide$last_deviation == 0, na.rm = TRUE)
if (.n_zero_dev > 0L) message("Figure S6: ", .n_zero_dev,
  " exact-zero deviations assigned to At or above zero.")
.n_nonfinite <- sum(!is.finite(sb_dev_wide$first_deviation) |
                      !is.finite(sb_dev_wide$last_deviation))
if (.n_nonfinite > 0L) {
  message("Figure S6: ", .n_nonfinite, " of ", nrow(sb_dev_wide),
          " RID x assay pairs have a non-finite deviation. Panel B excludes ",
          "them; panel A plots only the finite ones. The per-assay n printed ",
          "in panel A comes from the correlation step, which uses ",
          "complete observations.")
}
sb_write(sb_dev_wide, "Figure_S6_deviation_first_vs_last.csv")
sb_dev_rho_stat <- function(d) suppressWarnings(
  cor(d$first_deviation, d$last_deviation,
      method = "spearman", use = "complete.obs")
)
sb_dev_by_assay <- map_dfr(study_b_assays, function(a) {
  d <- filter(sb_dev_wide, assay == a)
  if (nrow(d) < 3L) {
    return(tibble(
      assay = a, n_participants = nrow(d),
      deviation_spearman = NA_real_,
      deviation_spearman_low = NA_real_,
      deviation_spearman_high = NA_real_
    ))
  }
  r <- sb_boot_summary(
    d, sb_dev_rho_stat, B = study_b_boot_n,
    seed = study_b_seed_dev + 10L
  )
  tibble(
    assay = a,
    n_participants = nrow(d),
    deviation_spearman = unname(r["estimate"]),
    deviation_spearman_low = unname(r["low"]),
    deviation_spearman_high = unname(r["high"])
  )
})
sb_write(sb_dev_by_assay, "Figure_S6_deviation_persistence_by_assay.csv")
assay_display <- function(x) {
  dplyr::recode(
    as.character(x),
    "AlzPath" = "ALZpath",
    "ALZpath" = "ALZpath",
    "C2N_ratio" = "C2N",
    "C2N" = "C2N",
    .default = as.character(x)
  )
}
FONT_FAMILY <- "sans"
COL_TEXT <- "#263238"
COL_AXIS <- "#65727A"
COL_POINT <- "#315F78"
COL_FIT <- "#D49A45"
COL_REFERENCE <- "#AAB2B7"
FIG_W_MM <- 190
FIG_H_MM <- 124
theme_fig3 <- sb_theme_classic(base_size = 8.5, base_family = FONT_FAMILY) +
  theme(
    text = element_text(colour = COL_TEXT),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.7, colour = COL_TEXT),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.line = element_line(colour = COL_AXIS, linewidth = 0.3),
    axis.ticks = element_line(colour = COL_AXIS, linewidth = 0.3),
    axis.ticks.length = grid::unit(1.2, "mm"),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    strip.background = element_blank(),
    strip.text = element_text(size = 8.2, face = "bold", colour = COL_TEXT,
                              margin = margin(b = 3)),
    plot.margin = margin(7, 8, 6, 7),
    plot.tag = element_text(size = 11, face = "bold", colour = COL_TEXT),
    plot.tag.position = "topleft",
    legend.position = "none"
  )
finite_dev <- c(sb_dev_wide$first_deviation, sb_dev_wide$last_deviation)
finite_dev <- finite_dev[is.finite(finite_dev)]
if (!length(finite_dev)) {
  stop("No finite first/last deviation pairs are available.", call. = FALSE)
}
dev_lim <- range(finite_dev)
dev_span <- diff(dev_lim)
if (!is.finite(dev_span) || dev_span == 0) dev_span <- 1
dev_lim <- dev_lim + c(-1, 1) * 0.04 * dev_span
sb_fig3a_labels <- sb_dev_by_assay %>%
  mutate(
    label = if_else(
      is.finite(deviation_spearman_low) & is.finite(deviation_spearman_high),
      sprintf("\u03c1 = %.2f [%.2f, %.2f]",
              deviation_spearman,
              deviation_spearman_low, deviation_spearman_high),
      sprintf("\u03c1 = %.2f", deviation_spearman)
    )
  )
p3a <- ggplot(sb_dev_wide, aes(first_deviation, last_deviation)) +
  geom_abline(slope = 1, intercept = 0, colour = COL_REFERENCE,
              linewidth = 0.45, linetype = "22") +
  geom_point(shape = 16, size = 1.25, alpha = 0.48, colour = COL_POINT) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = COL_FIT, linewidth = 0.65) +
  geom_text(
    data = sb_fig3a_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.08, vjust = 1.12, size = 2.45,
    lineheight = 0.95, family = FONT_FAMILY,
    colour = COL_TEXT, inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ assay, ncol = 2,
    labeller = ggplot2::as_labeller(assay_display)
  ) +
  coord_equal(xlim = dev_lim, ylim = dev_lim, expand = FALSE) +
  labs(
    x = "First-occasion adjusted deviation",
    y = "Last-occasion adjusted deviation",
    tag = "(A)"
  ) +
  theme_fig3 +
  theme(panel.spacing = grid::unit(5, "mm"))
sign_levels <- c("Below zero", "At or above zero")
sb_sign_transition <- sb_dev_wide %>%
  filter(is.finite(first_deviation), is.finite(last_deviation)) %>%
  transmute(
    assay = as.character(assay),
    first_sign = if_else(first_deviation < 0, "Below zero", "At or above zero"),
    last_sign = if_else(last_deviation < 0, "Below zero", "At or above zero")
  ) %>%
  count(assay, first_sign, last_sign, name = "n") %>%
  tidyr::complete(
    assay = as.character(study_b_assays),
    first_sign = sign_levels,
    last_sign = sign_levels,
    fill = list(n = 0L)
  ) %>%
  group_by(assay, first_sign) %>%
  mutate(
    first_total = sum(n),
    conditional_proportion = if_else(first_total > 0, n / first_total, NA_real_)
  ) %>%
  ungroup()
sb_sign_retention <- sb_sign_transition %>%
  group_by(assay) %>%
  summarise(
    n_participants = sum(n),
    retained_n = sum(n[first_sign == last_sign]),
    sign_retention = if_else(n_participants > 0,
                             retained_n / n_participants, NA_real_),
    .groups = "drop"
  )
.t4_path <- file.path(study_b_dir, "Table_S9_supporting_deviation_by_assay.csv")
.t4_soft <- identical(Sys.getenv("STUDY_B_FIG3_SKIP_TABLE4_CHECK", unset = ""), "1")
if (file.exists(.t4_path)) {
  .t4 <- sb_read_result(.t4_path, stringsAsFactors = FALSE)
  .t4_needed <- c("assay", "deviation_spearman", "sign_retention")
  .t4_missing <- setdiff(.t4_needed, names(.t4))
  if (length(.t4_missing)) {
    message("Figure S6: ", basename(.t4_path), " lacks column(s) ",
            paste(.t4_missing, collapse = ", "),
            "; the Table S9 supporting analysis cross-check was skipped.")
  } else {
    .chk <- sb_dev_by_assay %>%
      select(assay, fig_rho = deviation_spearman) %>%
      left_join(select(sb_sign_retention, assay, fig_ret = sign_retention),
                by = "assay") %>%
      left_join(
        .t4 %>% select(assay,
                       t4_rho = deviation_spearman,
                       t4_ret = sign_retention),
        by = "assay"
      ) %>%
      mutate(d_rho = abs(fig_rho - t4_rho),
             d_ret = abs(fig_ret - t4_ret))
    .unmatched <- filter(.chk, is.na(t4_rho) | is.na(t4_ret))
    .mismatch  <- filter(.chk,
                         (is.finite(d_rho) & d_rho > 1e-6) |
                           (is.finite(d_ret) & d_ret > 1e-6))
    .t4_age <- format(file.mtime(.t4_path), "%Y-%m-%d %H:%M:%S")
    if (nrow(.unmatched) || nrow(.mismatch)) {
      cat("\nFigure S6 versus Table S9 supporting analysis (", basename(.t4_path),
          ", written ", .t4_age, ")\n", sep = "")
      print(as.data.frame(.chk))
      .msg <- paste0(
        "Figure S6 does not reproduce Table S9 supporting analysis.\n",
        if (nrow(.unmatched))
          paste0("  Assay name(s) absent from Table S9 supporting analysis: ",
                 paste(.unmatched$assay, collapse = ", "), "\n") else "",
        if (nrow(.mismatch))
          paste0("  Numerical mismatch for: ",
                 paste(.mismatch$assay, collapse = ", "), "\n") else "",
        "Scripts 09 and 10 are running against different deviations. ",
        "Table S9 supporting analysis was written ", .t4_age, ". Re-run 03, then 11, then 12 in one ",
        "clean session, and do not use either version until they agree."
      )
      if (.t4_soft) {
        warning(.msg, call. = FALSE)
      } else {
        stop(.msg, call. = FALSE)
      }
    } else {
      message("Figure S6: per-assay Spearman rho and sign retention reproduce ",
              basename(.t4_path), " exactly (written ", .t4_age, ").")
    }
  }
} else {
  message("Figure S6: ", basename(.t4_path), " is not present, so the Table S9 supporting analysis ",
          "cross-check was skipped. Run script 10 in this session and re-run ",
          "this script to confirm the two agree.")
}
sb_sign_transition <- sb_sign_transition %>%
  left_join(sb_sign_retention, by = "assay") %>%
  mutate(
    transition = if_else(first_sign == last_sign,
                         "Retained sign", "Changed sign"),
    first_sign = factor(first_sign, levels = sign_levels),
    last_sign = factor(last_sign, levels = sign_levels),
    cell_label = if_else(
      is.finite(conditional_proportion),
      sprintf("%.0f%%\n(n = %d)", 100 * conditional_proportion, n),
      "NA"
    ),
    label_colour = if_else(
      is.finite(conditional_proportion) & conditional_proportion >= 0.55,
      "white", COL_TEXT
    )
  )
strip_lookup <- sb_sign_retention %>%
  mutate(
    assay_label = assay_display(assay),
    strip_label = if_else(
      is.finite(sign_retention),
      sprintf("%s\n%.1f%% retained", assay_label, 100 * sign_retention),
      sprintf("%s\nretention unavailable", assay_label)
    )
  )
strip_levels <- strip_lookup$strip_label[
  match(as.character(study_b_assays), strip_lookup$assay)
]
sb_sign_transition <- sb_sign_transition %>%
  left_join(select(strip_lookup, assay, strip_label), by = "assay") %>%
  mutate(strip_label = factor(strip_label, levels = strip_levels))
sb_write(
  sb_sign_transition %>%
    mutate(first_sign = as.character(first_sign),
           last_sign = as.character(last_sign),
           strip_label = as.character(strip_label)),
  "Figure_S6_deviation_sign_transition.csv"
)
p3b <- ggplot(
  sb_sign_transition,
  aes(last_sign, first_sign, fill = transition,
      alpha = conditional_proportion)
) +
  geom_tile(width = 0.92, height = 0.92,
            colour = "white", linewidth = 0.9) +
  geom_text(
    aes(label = cell_label, colour = label_colour),
    size = 2.25, lineheight = 0.92, family = FONT_FAMILY,
    alpha = 1
  ) +
  facet_wrap(~ strip_label, ncol = 2) +
  scale_fill_manual(
    values = c("Retained sign" = COL_POINT,
               "Changed sign" = COL_FIT),
    breaks = c("Retained sign", "Changed sign")
  ) +
  scale_alpha_continuous(range = c(0.18, 1), limits = c(0, 1),
                         guide = "none") +
  scale_colour_identity() +
  scale_x_discrete(labels = c("Below zero" = "Below",
                              "At or above zero" = "At/above"),
                   expand = expansion(add = 0.06)) +
  scale_y_discrete(labels = c("Below zero" = "Below",
                              "At or above zero" = "At/above"),
                   expand = expansion(add = 0.06)) +
  labs(
    x = "Last-occasion deviation sign",
    y = "First-occasion deviation sign",
    tag = "(B)"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                             override.aes = list(alpha = 1))) +
  theme_fig3 +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 7.2),
    strip.text = element_text(size = 7.8, lineheight = 0.95,
                              margin = margin(b = 3)),
    panel.spacing = grid::unit(4, "mm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.margin = margin(2, 0, 0, 0),
    legend.text = element_text(size = 7.3),
    legend.key.width = grid::unit(4.5, "mm"),
    legend.key.height = grid::unit(3, "mm"),
    plot.margin = margin(7, 7, 6, 9)
  )
fig3 <- patchwork::wrap_plots(p3a, p3b, nrow = 1, widths = c(1.48, 1.02)) +
  patchwork::plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 4, 3, 4)
    )
  )
out_fig3 <- file.path(study_b_dir, "Figure_S6_longitudinal_persistence")
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
sb_ggsave(paste0(out_fig3, ".pdf"), fig3,
       width = FIG_W_MM, height = FIG_H_MM, units = "mm",
       device = pdf_device, bg = "white")
sb_ggsave(paste0(out_fig3, ".tiff"), fig3,
       width = FIG_W_MM, height = FIG_H_MM, units = "mm",
       dpi = 600, compression = "lzw", bg = "white")
sb_ggsave(paste0(out_fig3, ".png"), fig3,
       width = FIG_W_MM, height = FIG_H_MM, units = "mm",
       dpi = 600, bg = "white")
sb_ggsave(file.path(study_b_dir, "Figure_S6_deviation_first_vs_last.png"), p3a,
       width = 125, height = 112, units = "mm", dpi = 600, bg = "white")
figS7_boot_ci <- function(d, seed) {
  person <- d %>%
    filter(!is.na(changed)) %>%
    group_by(RID) %>%
    summarise(value = mean(as.numeric(changed)), .groups = "drop")
  if (!nrow(person)) stop("Figure S7 has no complete participant data.", call. = FALSE)
  set.seed(seed)
  boot <- replicate(study_b_boot_n,
                    mean(sample(person$value, nrow(person), replace = TRUE)))
  tibble(
    comparisons = sum(!is.na(d$changed)),
    n_participants = nrow(person),
    estimate = mean(person$value),
    low = unname(quantile(boot, 0.025, na.rm = TRUE)),
    high = unname(quantile(boot, 0.975, na.rm = TRUE)),
    bootstrap_replicates = study_b_boot_n,
    estimand = "participant-equal mean; assay contrasts averaged within participant"
  )
}
same_long <- sb_dev_wide %>%
  filter(!is.na(first_class), !is.na(last_class)) %>%
  transmute(RID, contrast_id = as.character(assay),
            changed = first_class != last_class)
first_draw_wide <- sb_rep_dat %>%
  arrange(RID, assay, plasma_date, sample_id) %>%
  group_by(RID, assay) %>%
  summarise(first_class = first(as.character(three_class)), .groups = "drop") %>%
  pivot_wider(names_from = assay, values_from = first_class)
missing_assays <- setdiff(study_b_assays, names(first_draw_wide))
if (length(missing_assays)) {
  stop("Figure S7 is missing assay columns: ",
       paste(missing_assays, collapse = ", "), call. = FALSE)
}
assay_pairs <- combn(study_b_assays, 2, simplify = FALSE)
same_draw_switch <- map_dfr(assay_pairs, function(p) {
  first_draw_wide %>%
    transmute(RID, contrast_id = paste(p, collapse = " vs "),
              changed = .data[[p[1]]] != .data[[p[2]]])
})
sb_remeasure <- bind_rows(
  figS7_boot_ci(same_long, sb_seed("figS7_same_assay")) %>%
    mutate(contrast = "Same assay over time",
           display = "Same assay: first vs last draw"),
  figS7_boot_ci(same_draw_switch, sb_seed("figS7_same_draw_switch")) %>%
    mutate(contrast = "Different assays at same draw",
           display = "Different assays: same draw")
) %>%
  mutate(
    metric = "three_class_category_changed",
    per_100 = 100 * estimate,
    per_100_low = 100 * low,
    per_100_high = 100 * high
  ) %>%
  select(contrast, display, metric, comparisons, n_participants,
         per_100, per_100_low, per_100_high,
         bootstrap_replicates, estimand)
sb_write(sb_remeasure, "Figure_S7_remeasurement_vs_switching.csv")
dat_s5 <- sb_remeasure %>%
  mutate(
    display = factor(display, levels = rev(display)),
    value_label = sprintf("%.1f%% (95%% CI %.1f–%.1f)",
                          per_100, per_100_low, per_100_high)
  )
s5_xmax <- max(dat_s5$per_100_high, na.rm = TRUE)
s5_xmax <- max(10, ceiling((s5_xmax + 22) / 10) * 10)
pS5 <- ggplot(dat_s5, aes(per_100, display)) +
  geom_segment(aes(x = 0, xend = per_100,
                   y = display, yend = display),
               linewidth = 1.0, colour = "#C9D0D4", lineend = "round") +
  geom_errorbar(aes(xmin = per_100_low, xmax = per_100_high),
                width = 0.13, orientation = "y", colour = COL_POINT) +
  geom_point(shape = 21, size = 3.0, stroke = 0.8,
             fill = "white", colour = COL_POINT) +
  geom_text(aes(x = per_100_high + 1.2, label = value_label),
            hjust = 0, size = 2.65, family = FONT_FAMILY,
            colour = COL_TEXT) +
  scale_x_continuous(
    limits = c(0, s5_xmax),
    breaks = scales::pretty_breaks(5),
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "Three-class category changes per 100 comparisons", y = NULL) +
  theme_fig3 +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 8, margin = margin(r = 6)),
    plot.margin = margin(8, 32, 8, 8)
  )
out_s7 <- file.path(study_b_dir, "Figure_S7_remeasurement_vs_switching")
sb_ggsave(paste0(out_s7, ".pdf"), pS5,
       width = 165, height = 74, units = "mm", device = pdf_device, bg = "white")
sb_ggsave(paste0(out_s7, ".tiff"), pS5,
       width = 165, height = 74, units = "mm", dpi = 600,
       compression = "lzw", bg = "white")
sb_ggsave(paste0(out_s7, ".png"), pS5,
       width = 165, height = 74, units = "mm", dpi = 600, bg = "white")
cat("\nFigure S6 -- per-assay deviation persistence\n")
print(sb_dev_by_assay)
cat("\nFigure S6 -- per-assay sign retention\n")
print(sb_sign_retention)
cat("\nFigure S7 -- like-for-like three-class comparison\n")
print(sb_remeasure)
