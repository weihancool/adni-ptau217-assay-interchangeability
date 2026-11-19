if (!exists("sb_begin_reporting_step", mode = "function"))
  stop("Run revised 01_config.R first.", call. = FALSE)
sb_begin_reporting_step("16_figure2_switching_consequences.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
if (exists("study_b_dir", inherits = TRUE)) {
  res <- get("study_b_dir", inherits = TRUE)
} else {
  root <- Sys.getenv("STUDY_B_ROOT")
  if (!nzchar(root)) {
    stop(
      "Cannot locate RESULTS. Set STUDY_B_ROOT to your project folder, ",
      "or use the R session where study_b_dir is already defined.",
      call. = FALSE
    )
  }
  res <- file.path(root, "RESULTS")
}
if (!is.character(res) || length(res) != 1L ||
    is.na(res) || !dir.exists(res)) {
  stop("RESULTS folder does not exist: ", res, call. = FALSE)
}
fig2_read_frozen <- function(path) {
  sb_read_result(path, check.names = FALSE, stringsAsFactors = FALSE)
}
a1 <- fig2_read_frozen(
  file.path(res, "Table_S7_A1_assay_switching_consequences.csv")
)
margin <- fig2_read_frozen(
  file.path(res, "Table_S7_A1_switching_risk_by_source_margin_CI.csv")
)
a2 <- fig2_read_frozen(
  file.path(res, "Table_S8_A2_symptomatic_switching_consequences.csv")
)
require_columns <- function(z, columns, context) {
  missing_columns <- setdiff(columns, names(z))
  if (length(missing_columns)) {
    stop(
      context, ": missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}
require_columns(
  a1,
  c(
    "population", "scope", "metric",
    "from_assay", "to_assay",
    "per_100", "per_100_low", "per_100_high"
  ),
  "A1 results"
)
require_columns(
  margin,
  c(
    "metric", "source_margin_bin",
    "per_100", "per_100_low", "per_100_high"
  ),
  "Source-margin results"
)
require_columns(
  a2,
  c(
    "metric", "from_assay", "to_assay",
    "per_100", "per_100_low", "per_100_high"
  ),
  "A2 results"
)
check_rates <- function(z, context) {
  columns <- c("per_100", "per_100_low", "per_100_high")
  require_columns(z, columns, context)
  if (!nrow(z) ||
      !all(vapply(z[columns], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(z[columns])))) {
    stop(
      "Figure 2: missing or invalid estimates for ", context,
      call. = FALSE
    )
  }
  if (any(z$per_100_low > z$per_100_high) ||
      any(as.matrix(z[columns]) < 0) ||
      any(as.matrix(z[columns]) > 100)) {
    stop(
      "Figure 2: invalid rates or confidence intervals for ",
      context,
      call. = FALSE
    )
  }
}
assays <- c("AlzPath", "C2N_ratio", "Fujirebio", "Janssen")
assay_labels <- c("ALZpath", "C2N", "Fujirebio", "Janssen")
ink <- "#253444"
blue <- "#254F70"
red <- "#B44145"
grid_col <- "#E7EBEE"
base_theme <- theme_classic(base_size = 9) +
  theme(
    text = element_text(colour = ink),
    axis.text = element_text(colour = ink),
    axis.title = element_text(colour = ink),
    axis.line = element_line(colour = ink, linewidth = 0.35),
    axis.ticks = element_line(colour = ink, linewidth = 0.35),
    plot.title = element_text(
      face = "bold", size = 12, hjust = 0
    ),
    legend.key = element_blank(),
    plot.margin = margin(5.5, 7, 5.5, 5.5)
  )
pair <- a1 %>%
  filter(
    population == "Full cohort",
    scope == "Primary anchor",
    metric == "binary_flip",
    from_assay != "All",
    to_assay != "All"
  )
if (nrow(pair) != 12L ||
    anyNA(pair[c("from_assay", "to_assay")]) ||
    anyDuplicated(pair[c("from_assay", "to_assay")]) ||
    any(pair$from_assay == pair$to_assay) ||
    any(!pair$from_assay %in% assays) ||
    any(!pair$to_assay %in% assays)) {
  stop(
    "Figure 2 requires 12 unique off-diagonal assay pairs.",
    call. = FALSE
  )
}
check_rates(pair, "panel A")
pair <- pair %>%
  mutate(
    from_assay = factor(from_assay, levels = rev(assays)),
    to_assay = factor(to_assay, levels = assays)
  )
pA <- ggplot(pair, aes(to_assay, from_assay, fill = per_100)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf("%.1f", per_100)),
    size = 2.8,
    colour = ink
  ) +
  scale_x_discrete(
    limits = assays,
    labels = assay_labels,
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(assays),
    labels = rev(assay_labels),
    drop = FALSE
  ) +
  scale_fill_gradient(
    low = "#EDF1F4",
    high = "#C46565",
    name = "Reclassifications per 100"
  ) +
  coord_equal() +
  labs(
    title = "A",
    x = "Destination assay",
    y = "Source assay"
  ) +
  base_theme +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7)
  ) +
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(26, "mm"),
      barheight = grid::unit(2, "mm"),
      title.position = "top"
    )
  )
mr <- margin %>%
  filter(metric == "binary_flip")
check_rates(mr, "panel B")
if (anyNA(mr$source_margin_bin) ||
    any(!nzchar(as.character(mr$source_margin_bin))) ||
    anyDuplicated(mr$source_margin_bin)) {
  stop(
    "Figure 2: missing or duplicate source-margin bins.",
    call. = FALSE
  )
}
bin_order <- as.character(mr$source_margin_bin)
if ("participants" %in% names(mr)) {
  if (!is.numeric(mr$participants) ||
      any(!is.finite(mr$participants)) ||
      any(mr$participants < 1) ||
      any(mr$participants != floor(mr$participants))) {
    stop(
      "Figure 2: invalid participant counts in panel B.",
      call. = FALSE
    )
  }
  bin_labels <- paste0(
    bin_order,
    "  (n = ",
    mr$participants,
    ")"
  )
} else {
  bin_labels <- bin_order
}
mr$bin_display <- factor(
  bin_order,
  levels = rev(bin_order)
)
metric_labels <- c(
  binary_flip = "Binary reclassification",
  three_class_change = "Three-class change",
  confirmatory_change = "Change involving\nintermediate category"
)
overall <- bind_rows(
  a1 %>%
    filter(
      population == "Full cohort",
      scope == "Primary anchor",
      from_assay == "All",
      to_assay == "All"
    ) %>%
    mutate(group = "Full cohort"),
  a2 %>%
    filter(
      from_assay == "All",
      to_assay == "All"
    ) %>%
    mutate(group = "MCI or dementia")
) %>%
  filter(metric %in% names(metric_labels)) %>%
  mutate(
    group = factor(
      group,
      levels = c("Full cohort", "MCI or dementia")
    ),
    metric_label = factor(
      unname(metric_labels[metric]),
      levels = rev(unname(metric_labels))
    )
  )
check_rates(overall, "panel C")
if (nrow(overall) != 6L ||
    anyDuplicated(overall[c("group", "metric")])) {
  stop(
    "Figure 2 requires three metrics per population.",
    call. = FALSE
  )
}
rate_max <- min(
  100,
  max(
    10,
    ceiling(
      max(c(mr$per_100_high, overall$per_100_high)) / 10
    ) * 10
  )
)
rate_breaks <- pretty(c(0, rate_max), n = 4)
rate_breaks <- rate_breaks[
  rate_breaks >= 0 & rate_breaks <= rate_max
]
pB <- ggplot(mr, aes(per_100, bin_display)) +
  geom_errorbar(
    aes(xmin = per_100_low, xmax = per_100_high),
    orientation = "y",
    width = 0.16,
    linewidth = 0.5,
    colour = blue
  ) +
  geom_point(size = 2, colour = blue) +
  scale_y_discrete(
    labels = setNames(bin_labels, bin_order)
  ) +
  scale_x_continuous(
    limits = c(0, rate_max),
    breaks = rate_breaks,
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "B",
    x = "Binary reclassifications per 100 switches",
    y = "Absolute source margin (SD units)"
  ) +
  base_theme +
  theme(
    panel.grid.major.x = element_line(
      colour = grid_col,
      linewidth = 0.3
    ),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(size = 7)
  )
pd <- position_dodge(width = 0.38, orientation = "y")
pC <- ggplot(
  overall,
  aes(
    per_100,
    metric_label,
    colour = group,
    shape = group,
    group = group
  )
) +
  geom_errorbar(
    aes(xmin = per_100_low, xmax = per_100_high),
    orientation = "y",
    width = 0.13,
    linewidth = 0.6,
    position = pd
  ) +
  geom_point(size = 2.3, position = pd) +
  scale_colour_manual(
    values = c(
      "Full cohort" = blue,
      "MCI or dementia" = red
    )
  ) +
  scale_shape_manual(
    values = c(
      "Full cohort" = 16,
      "MCI or dementia" = 17
    )
  ) +
  scale_x_continuous(
    limits = c(0, rate_max),
    breaks = rate_breaks,
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "C",
    x = "Events per 100 directed switches",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  base_theme +
  theme(
    panel.grid.major.x = element_line(
      colour = grid_col,
      linewidth = 0.3
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8)
  )
fig <- ((pA | pB) / pC) +
  plot_layout(heights = c(1.35, 0.8))
base <- file.path(
  res,
  "Figure_2_switching_consequences"
)
pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}
sb_ggsave(
  filename = paste0(base, ".pdf"),
  plot = fig,
  device = pdf_device,
  width = 190,
  height = 160,
  units = "mm",
  bg = "white"
)
sb_ggsave(
  filename = paste0(base, ".png"),
  plot = fig,
  device = "png",
  width = 190,
  height = 160,
  units = "mm",
  dpi = 600,
  bg = "white"
)
sb_ggsave(
  filename = paste0(base, ".tiff"),
  plot = fig,
  device = "tiff",
  width = 190,
  height = 160,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
source_data <- bind_rows(
  pair %>%
    mutate(
      from_assay = as.character(from_assay),
      to_assay = as.character(to_assay),
      panel = "A"
    ),
  mr %>%
    mutate(
      bin_display = as.character(bin_display),
      panel = "B"
    ),
  overall %>%
    mutate(
      group = as.character(group),
      metric_label = as.character(metric_label),
      panel = "C"
    )
)
sb_save_csv(
  source_data,
  file = file.path(res, "Figure_2_source_data.csv"),
  row.names = FALSE,
  na = ""
)
caption <- paste0(
  "Figure 2. Same-draw counterfactual assay switching. ",
  "A: binary reclassification per 100 switches for each ",
  "source-to-destination pair in primary anchors; diagonal ",
  "cells are not switches. ",
  "B: pooled binary reclassification by absolute source margin ",
  "across all collection occasions. Points and bars show ",
  "estimates and participant-cluster bootstrap 95% CIs. ",
  "n denotes unique participants in each bin, not independent ",
  "switches; participants can occur in multiple bins. ",
  "C: consequences in primary anchors of the full cohort and ",
  "MCI/dementia subgroup; points and bars show equal-weight ",
  "averages across 12 directed switches and their existing ",
  "95% CIs. ",
  "Change involving intermediate category denotes a change ",
  "between the intermediate category and either the negative ",
  "or positive category, excluding direct negative-to-positive ",
  "or positive-to-negative switches. ",
  "The binary pairwise event is symmetric; the 12 directed cells correspond to six unordered assay pairs. ",
  "These are counterfactual switches, not observed changes ",
  "in clinical practice."
)
writeLines(
  caption,
  con = file.path(res, "Figure_2_caption.txt"),
  useBytes = TRUE
)
message(
  "Figure 2 saved in: ", normalizePath(res),
  "\nOutputs: Figure_2_switching_consequences.pdf / png / tiff",
  "\nExisting estimates and confidence intervals were retained."
)
sb_stamp_result(file.path(res, "Figure_2_caption.txt"))
sb_end_reporting_step()
if (interactive()) {
  print(fig)
}
