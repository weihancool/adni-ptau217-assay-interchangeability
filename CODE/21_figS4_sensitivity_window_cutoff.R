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
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "purrr", "lme4", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Please install missing package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
})
if (!exists("sb_build_events", inherits = TRUE)) {
  stop("Run 02_preprocess.R manually first.", call. = FALSE)
}
if (!exists("sb_crossfit", inherits = TRUE)) {
  stop("Run 03_core_computation.R manually first.", call. = FALSE)
}
required_objects <- c(
  "study_b_sens_windows", "study_b_sens_cutoffs", "study_b_assays",
  "amy_pre", "sb_checked_lmer", "sb_safe_ratio", "sb_write"
)
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects)) {
  stop(
    "Required Study B object(s) unavailable: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}
if (!all(c("RID", "date", "qc_flag_raw", "amyloid_centiloid") %in%
         names(amy_pre))) {
  stop(
    "amy_pre must contain RID, date, qc_flag_raw, and amyloid_centiloid.",
    call. = FALSE
  )
}
if (length(study_b_sens_windows) < 1L ||
    length(study_b_sens_cutoffs) < 1L) {
  stop("Sensitivity windows and cutoffs must be non-empty.", call. = FALSE)
}
assay_order <- as.character(study_b_assays)
if (!length(assay_order) || anyDuplicated(assay_order)) {
  stop("study_b_assays must contain unique assay names.", call. = FALSE)
}
sb_sens_grid <- tidyr::expand_grid(
  window_days = sort(unique(as.numeric(study_b_sens_windows))),
  centiloid_cutoff = sort(unique(as.numeric(study_b_sens_cutoffs))),
  pet_qc_rule = c("full_pass", "full_or_partial")
)
if (any(!is.finite(sb_sens_grid$window_days)) ||
    any(!is.finite(sb_sens_grid$centiloid_cutoff))) {
  stop("Sensitivity windows and cutoffs must be finite numeric values.",
       call. = FALSE)
}
.sb_sens_status_row <- function(
    window_days, centiloid_cutoff, pet_qc_rule, status) {
  tibble(
    window_days = as.numeric(window_days),
    centiloid_cutoff = as.numeric(centiloid_cutoff),
    pet_qc_rule = as.character(pet_qc_rule),
    status = as.character(status),
    matched_events = NA_integer_,
    matched_participants = NA_integer_,
    anchor_n = NA_integer_,
    anchor_pet_available_n = NA_integer_,
    amyloid_positive_anchor_n = NA_integer_,
    youden_AARR = NA_real_,
    youden_AARR_complete_n = NA_integer_,
    three_class_disagreement = NA_real_,
    three_class_complete_n = NA_integer_,
    person_by_assay_proportion = NA_real_,
    occasion_proportion = NA_real_,
    assay_as_fixed = TRUE,
    person_assay_to_occasion_ratio = NA_real_,
    ratio_definition =
      "person_by_assay / shared_occasion_within_person",
    singular = NA
  )
}
.sb_discordant <- function(m) {
  m <- m[rowSums(is.na(m)) == 0L, , drop = FALSE]
  if (!nrow(m)) return(c(rate = NA_real_, n = 0))
  c(
    rate = mean(apply(m, 1, function(z) length(unique(z)) > 1L)),
    n = nrow(m)
  )
}
sb_sens_one <- function(window_days, centiloid_cutoff, pet_qc_rule) {
  message(
    "Figure S4: window=", window_days,
    " days; cutoff=", centiloid_cutoff,
    "; PET QC=", pet_qc_rule
  )
  pet_levels <- if (identical(pet_qc_rule, "full_pass")) {
    2L
  } else {
    c(1L, 2L)
  }
  eligible_pet <- amy_pre %>%
    filter(
      qc_flag_raw %in% pet_levels,
      is.finite(amyloid_centiloid)
    )
  conflicts <- eligible_pet %>%
    group_by(RID, date) %>%
    summarise(
      values = n_distinct(amyloid_centiloid),
      .groups = "drop"
    ) %>%
    filter(values > 1L)
  if (nrow(conflicts)) {
    return(.sb_sens_status_row(
      window_days, centiloid_cutoff, pet_qc_rule,
      paste0("failed: conflicting eligible PET values in ",
             nrow(conflicts), " RID-date cells")
    ))
  }
  pet_sens <- eligible_pet %>%
    arrange(RID, date, desc(qc_flag_raw)) %>%
    distinct(RID, date, .keep_all = TRUE)
  ev <- tryCatch(
    sb_build_events(
      window_days = window_days,
      centiloid_cutoff = centiloid_cutoff,
      pet_data = pet_sens
    ),
    error = function(e) e
  )
  if (inherits(ev, "error")) {
    return(.sb_sens_status_row(
      window_days, centiloid_cutoff, pet_qc_rule,
      paste("failed: event construction:", conditionMessage(ev))
    ))
  }
  required_event_columns <- c("RID", "amyloid_positive")
  if (!all(required_event_columns %in% names(ev))) {
    return(.sb_sens_status_row(
      window_days, centiloid_cutoff, pet_qc_rule,
      "failed: event data missing RID or amyloid_positive"
    ))
  }
  ok <- ev %>% filter(!is.na(amyloid_positive))
  if (n_distinct(ok$RID) < 100L ||
      n_distinct(ok$amyloid_positive) < 2L) {
    return(.sb_sens_status_row(
      window_days, centiloid_cutoff, pet_qc_rule,
      "insufficient reference data"
    ))
  }
  sc <- tryCatch(sb_crossfit(ev), error = function(e) e)
  if (inherits(sc, "error")) {
    return(.sb_sens_status_row(
      window_days, centiloid_cutoff, pet_qc_rule,
      paste("failed: cross-fitting:", conditionMessage(sc))
    ))
  }
  required_sc_columns <- c(
    "RID", "sample_id", "assay", "primary_anchor", "amyloid_positive",
    "youden_positive", "three_class", "margin_youden",
    "years_from_baseline", "baseline_amyloid_group", "phase"
  )
  missing_sc_columns <- setdiff(required_sc_columns, names(sc))
  if (length(missing_sc_columns)) {
    return(.sb_sens_status_row(
      window_days, centiloid_cutoff, pet_qc_rule,
      paste(
        "failed: cross-fit output missing",
        paste(missing_sc_columns, collapse = ", ")
      )
    ))
  }
  anch <- sc %>% filter(primary_anchor %in% TRUE)
  if (!nrow(anch)) {
    return(.sb_sens_status_row(
      window_days, centiloid_cutoff, pet_qc_rule,
      "failed: no primary-anchor observations"
    ))
  }
  duplicated_anchor_cells <- anch %>%
    count(RID, assay, name = "records") %>%
    filter(records > 1L)
  if (nrow(duplicated_anchor_cells)) {
    return(.sb_sens_status_row(
      window_days, centiloid_cutoff, pet_qc_rule,
      paste0(
        "failed: duplicated primary-anchor RID-assay cells (",
        nrow(duplicated_anchor_cells), ")"
      )
    ))
  }
  w <- anch %>%
    select(RID, assay, value = youden_positive) %>%
    pivot_wider(names_from = assay, values_from = value)
  wc <- anch %>%
    select(RID, assay, three_class) %>%
    mutate(three_class = as.character(three_class)) %>%
    pivot_wider(names_from = assay, values_from = three_class)
  missing_binary_assays <- setdiff(assay_order, names(w))
  missing_class_assays <- setdiff(assay_order, names(wc))
  if (length(missing_binary_assays) || length(missing_class_assays)) {
    return(.sb_sens_status_row(
      window_days, centiloid_cutoff, pet_qc_rule,
      paste0(
        "failed: missing assay columns; binary=",
        paste(missing_binary_assays, collapse = "/"),
        "; three-class=",
        paste(missing_class_assays, collapse = "/")
      )
    ))
  }
  aarr_r <- .sb_discordant(
    as.matrix(w[, assay_order, drop = FALSE])
  )
  three_r <- .sb_discordant(
    as.matrix(wc[, assay_order, drop = FALSE])
  )
  d <- sc %>%
    mutate(
      occ_id = interaction(RID, sample_id, drop = TRUE),
      phase = factor(phase),
      assay = factor(assay, levels = assay_order),
      RID = factor(RID)
    )
  model_fit <- tryCatch(
    sb_checked_lmer(
      as.formula(paste("margin_youden ~",
        paste(c(sb_vc_fixed_terms(d), "(1|RID)", "(1|occ_id)", "(1|RID:assay)"),
              collapse = " + "))),
      data = d,
      REML = TRUE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) e
  )
  if (inherits(model_fit, "error")) {
    return(.sb_sens_status_row(
      window_days, centiloid_cutoff, pet_qc_rule,
      paste("failed: variance model:", conditionMessage(model_fit))
    ))
  }
  v <- as.data.frame(lme4::VarCorr(model_fit))
  pa <- v$vcov[v$grp %in% c("RID:assay", "assay:RID")][1]
  po <- v$vcov[v$grp == "occ_id"][1]
  tot <- sum(v$vcov[is.finite(v$vcov)])
  if (!is.finite(pa) || !is.finite(po) || !is.finite(tot) || tot <= 0) {
    return(.sb_sens_status_row(
      window_days, centiloid_cutoff, pet_qc_rule,
      "failed: required variance components unavailable"
    ))
  }
  tibble(
    window_days = as.numeric(window_days),
    centiloid_cutoff = as.numeric(centiloid_cutoff),
    pet_qc_rule = as.character(pet_qc_rule),
    status = "ok",
    matched_events = n_distinct(
      sc$sample_id[!is.na(sc$amyloid_positive)]
    ),
    matched_participants = n_distinct(
      sc$RID[!is.na(sc$amyloid_positive)]
    ),
    anchor_n = n_distinct(anch$RID),
    anchor_pet_available_n = n_distinct(
      anch$RID[!is.na(anch$amyloid_positive)]
    ),
    amyloid_positive_anchor_n = n_distinct(
      anch$RID[
        anch$amyloid_positive == 1 & !is.na(anch$amyloid_positive)
      ]
    ),
    youden_AARR = unname(aarr_r["rate"]),
    youden_AARR_complete_n = as.integer(unname(aarr_r["n"])),
    three_class_disagreement = unname(three_r["rate"]),
    three_class_complete_n = as.integer(unname(three_r["n"])),
    person_by_assay_proportion = pa / tot,
    occasion_proportion = po / tot,
    assay_as_fixed = TRUE,
    person_assay_to_occasion_ratio = sb_safe_ratio(pa, po),
    ratio_definition =
      "person_by_assay / shared_occasion_within_person",
    singular = lme4::isSingular(model_fit)
  )
}
sb_sensitivity <- purrr::pmap_dfr(sb_sens_grid, sb_sens_one)
sb_write(sb_sensitivity, "Figure_S4_sensitivity_window_cutoff.csv")
sb_write(sb_sensitivity, "Table_S6_sensitivity_window_cutoff_full.csv")
if (identical(as.integer(study_b_pet_qc_allowed), 2L) &&
    sb_have("sb_scored", "sb_anchor")) {
  .sb13_primary <- sb_sensitivity %>%
    filter(window_days == study_b_amyloid_window,
           centiloid_cutoff == study_b_centiloid_cutoff,
           pet_qc_rule == "full_pass")
  .sb13_expected <- c(
    matched_events = n_distinct(sb_scored$sample_id[!is.na(sb_scored$amyloid_positive)]),
    matched_participants = n_distinct(sb_scored$RID[!is.na(sb_scored$amyloid_positive)]),
    amyloid_positive_anchor_n = n_distinct(sb_anchor$RID[sb_anchor$amyloid_positive %in% 1L]),
    anchor_pet_available_n = n_distinct(sb_anchor$RID[!is.na(sb_anchor$amyloid_positive)]),
    youden_AARR = sb_discordant_rate(sb_rule_matrix(sb_anchor, "youden_positive")))
  if (nrow(.sb13_primary) != 1L || !identical(.sb13_primary$status, "ok")) {
    stop("Figure S4: primary specification row missing or failed.", call. = FALSE)
  }
  .sb13_observed <- unlist(.sb13_primary[1, names(.sb13_expected)])
  .sb13_bad <- names(.sb13_expected)[!(abs(.sb13_observed - .sb13_expected) <= 1e-8)]
  if (length(.sb13_bad)) {
    stop("Figure S4: primary sensitivity row does not reproduce the main analysis (",
         paste(.sb13_bad, collapse = ", "), "). Rerun 01-13 in one clean session.",
         call. = FALSE)
  }
  message("Figure S4: primary specification reproduces the main analysis.")
}
sb_sensitivity_ok <- sb_sensitivity %>%
  filter(status == "ok")
if (nrow(sb_sensitivity_ok)) {
  window_levels <- sort(unique(sb_sensitivity_ok$window_days))
  window_labels <- paste0(window_levels, " days")
  sb_plot <- sb_sensitivity_ok %>%
    select(
      window_days, centiloid_cutoff, pet_qc_rule,
      person_assay_to_occasion_ratio,
      youden_AARR,
      three_class_disagreement
    ) %>%
    pivot_longer(
      cols = c(
        person_assay_to_occasion_ratio,
        youden_AARR,
        three_class_disagreement
      ),
      names_to = "metric",
      values_to = "value"
    ) %>%
    filter(is.finite(value)) %>%
    mutate(
      metric_label = recode(
        metric,
        "person_assay_to_occasion_ratio" = "Variance ratio",
        "youden_AARR" = "Youden disagreement (%)",
        "three_class_disagreement" = "Three-class disagreement (%)"
      ),
      metric_label = factor(
        metric_label,
        levels = c(
          "Variance ratio",
          "Youden disagreement (%)",
          "Three-class disagreement (%)"
        )
      ),
      value_plot = if_else(
        metric == "person_assay_to_occasion_ratio",
        value,
        100 * value
      ),
      qc_label = recode(
        pet_qc_rule,
        "full_pass" = "Full-pass visual QC",
        "full_or_partial" = "Full + partial-pass visual QC"
      ),
      qc_label = factor(
        qc_label,
        levels = c(
          "Full-pass visual QC",
          "Full + partial-pass visual QC"
        )
      ),
      window_label = factor(
        paste0(window_days, " days"),
        levels = window_labels
      )
    )
  ratio_reference <- tibble(
    metric_label = factor(
      "Variance ratio",
      levels = levels(sb_plot$metric_label)
    ),
    reference = 1
  )
  FONT_FAMILY <- "sans"
  COL_TEXT <- "#263238"
  COL_AXIS <- "#65727A"
  COL_GRID <- "#E7ECEF"
  COL_REFERENCE <- "#9CA7AD"
  WINDOW_COLOURS <- c("#315F78", "#D49A45", "#7B6FA6", "#4F8A73")
  WINDOW_SHAPES <- c(21, 24, 22, 23)
  WINDOW_LINES <- c("solid", "32", "42", "dotdash")
  FIG_W_MM <- 183
  FIG_H_MM <- 148
  n_windows <- length(window_labels)
  if (n_windows > length(WINDOW_COLOURS)) {
    stop(
      "Figure S4 supports at most ", length(WINDOW_COLOURS),
      " PET-window values without extending the plotting palette.",
      call. = FALSE
    )
  }
  window_colours <- setNames(WINDOW_COLOURS[seq_len(n_windows)],
                             window_labels)
  window_shapes <- setNames(WINDOW_SHAPES[seq_len(n_windows)],
                            window_labels)
  window_lines <- setNames(WINDOW_LINES[seq_len(n_windows)],
                           window_labels)
  position_windows <- position_dodge(width = 0.65)
  p <- ggplot(
    sb_plot,
    aes(
      x = centiloid_cutoff,
      y = value_plot,
      group = window_label,
      colour = window_label,
      shape = window_label,
      linetype = window_label
    )
  ) +
    geom_hline(
      data = ratio_reference,
      aes(yintercept = reference),
      inherit.aes = FALSE,
      colour = COL_REFERENCE,
      linewidth = 0.42,
      linetype = "22"
    ) +
    geom_line(
      linewidth = 0.62,
      position = position_windows
    ) +
    geom_point(
      size = 2.35,
      stroke = 0.72,
      fill = "white",
      position = position_windows
    ) +
    facet_grid(
      rows = vars(metric_label),
      cols = vars(qc_label),
      scales = "free_y"
    ) +
    scale_x_continuous(
      breaks = sort(unique(sb_plot$centiloid_cutoff)),
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 4),
      expand = expansion(mult = c(0.08, 0.10))
    ) +
    scale_colour_manual(
      values = window_colours,
      breaks = window_labels,
      name = "PET window"
    ) +
    scale_shape_manual(
      values = window_shapes,
      breaks = window_labels,
      name = "PET window"
    ) +
    scale_linetype_manual(
      values = window_lines,
      breaks = window_labels,
      name = "PET window"
    ) +
    labs(
      x = "Centiloid positivity cutoff",
      y = NULL,
      caption = paste0(
        "Window series are offset slightly along the x-axis for visibility; ",
        "equal vertical positions indicate identical estimates. ",
        "Each participant contributes the nearest eligible PET reference."
      )
    ) +
    guides(
      colour = guide_legend(
        order = 1, nrow = 1,
        override.aes = list(fill = "white", linewidth = 0.7)
      ),
      shape = guide_legend(order = 1, nrow = 1),
      linetype = guide_legend(order = 1, nrow = 1)
    ) +
    sb_theme_classic(base_size = 8, base_family = FONT_FAMILY) +
    theme(
      text = element_text(colour = COL_TEXT),
      panel.grid.major.y = element_line(
        colour = COL_GRID,
        linewidth = 0.28
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = COL_AXIS, linewidth = 0.32),
      axis.ticks = element_line(colour = COL_AXIS, linewidth = 0.32),
      axis.ticks.length = grid::unit(1.1, "mm"),
      axis.text = element_text(size = 7.2, colour = COL_TEXT),
      axis.title.x = element_text(size = 8, margin = margin(t = 7)),
      strip.background = element_blank(),
      strip.text.x = element_text(
        size = 8, face = "bold", colour = COL_TEXT,
        margin = margin(b = 4)
      ),
      strip.text.y.right = element_text(
        size = 7.2, colour = COL_TEXT,
        angle = 0, hjust = 0,
        margin = margin(l = 5, r = 2)
      ),
      panel.spacing.x = grid::unit(7, "mm"),
      panel.spacing.y = grid::unit(5, "mm"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.title = element_text(size = 7.3),
      legend.text = element_text(size = 7.2),
      legend.key.width = grid::unit(6, "mm"),
      legend.key.height = grid::unit(3.4, "mm"),
      legend.margin = margin(t = 2, unit = "pt"),
      plot.caption = element_text(
        size = 6.6, colour = COL_AXIS,
        hjust = 0, lineheight = 1.05,
        margin = margin(t = 6)
      ),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(8, 8, 7, 8)
    )
  output_base <- file.path(study_b_dir, "Figure_S4_sensitivity_ratio")
  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }
  sb_ggsave(
    paste0(output_base, ".pdf"), p,
    width = FIG_W_MM, height = FIG_H_MM, units = "mm",
    device = pdf_device, bg = "white"
  )
  sb_ggsave(
    paste0(output_base, ".tiff"), p,
    width = FIG_W_MM, height = FIG_H_MM, units = "mm",
    dpi = 600, compression = "lzw", bg = "white"
  )
  sb_ggsave(
    paste0(output_base, ".png"), p,
    width = FIG_W_MM, height = FIG_H_MM, units = "mm",
    dpi = 600, bg = "white"
  )
} else {
  message(
    "Figure S4: every specification failed; plot skipped. ",
    "Inspect status in Figure_S4_sensitivity_window_cutoff.csv."
  )
}
cat("\nFigure S4 -- window × cutoff × PET-QC sensitivity\n")
print(as.data.frame(sb_sensitivity))
