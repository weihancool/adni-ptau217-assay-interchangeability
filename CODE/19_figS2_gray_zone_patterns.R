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
required_packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "scales")
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
  library(patchwork)
})
if (!exists("sb_anchor", inherits = TRUE)) {
  stop("Run 03_core_computation.R manually first.", call. = FALSE)
}
required_columns <- c("RID", "assay", "three_class")
missing_columns <- setdiff(required_columns, names(sb_anchor))
if (length(missing_columns)) {
  stop(
    "sb_anchor is missing: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}
if (!exists("study_b_assays", inherits = TRUE) || !length(study_b_assays)) {
  stop("study_b_assays is unavailable or empty.", call. = FALSE)
}
assay_order <- as.character(study_b_assays)
if (anyDuplicated(assay_order)) {
  stop("study_b_assays contains duplicated assay names.", call. = FALSE)
}
if (length(assay_order) != 4L) {
  stop("Figure S2 requires exactly four assays in study_b_assays.",
       call. = FALSE)
}
canonical_three_class <- function(x) {
  raw <- trimws(tolower(as.character(x)))
  token <- gsub("[^a-z0-9]+", " ", raw)
  token <- trimws(gsub("[[:space:]]+", " ", token))
  out <- rep(NA_character_, length(token))
  out[token %in% c("0", "negative", "neg", "below", "low", "normal")] <-
    "Negative"
  out[token %in% c(
    "1", "intermediate", "gray", "grey", "gray zone", "grey zone",
    "indeterminate", "borderline", "middle"
  )] <- "Intermediate"
  out[token %in% c("2", "positive", "pos", "above", "high", "abnormal")] <-
    "Positive"
  unresolved <- is.na(out) & !is.na(raw) & nzchar(raw)
  out[unresolved & grepl(
    "(^| )(gray|grey|intermediate|indeterminate|borderline|adjacent)( |$)",
    token
  )] <- "Intermediate"
  unresolved <- is.na(out) & !is.na(raw) & nzchar(raw)
  out[unresolved & grepl(
    "(^| )(negative|neg|below|low)( |$)", token
  )] <- "Negative"
  unresolved <- is.na(out) & !is.na(raw) & nzchar(raw)
  out[unresolved & grepl(
    "(^| )(positive|pos|above|high)( |$)", token
  )] <- "Positive"
  out
}
sb_class_long <- sb_anchor %>%
  transmute(
    RID,
    assay = as.character(assay),
    class_raw = as.character(three_class),
    class = canonical_three_class(three_class)
  ) %>%
  filter(assay %in% assay_order)
unmapped_classes <- sb_class_long %>%
  filter(!is.na(class_raw), nzchar(trimws(class_raw)), is.na(class)) %>%
  distinct(class_raw) %>%
  pull(class_raw)
if (length(unmapped_classes)) {
  stop(
    "Unrecognized three_class label(s): ",
    paste(unmapped_classes, collapse = ", "),
    ". Update canonical_three_class() before plotting.",
    call. = FALSE
  )
}
duplicate_cells <- sb_class_long %>%
  count(RID, assay, name = "records") %>%
  filter(records > 1L)
if (nrow(duplicate_cells)) {
  stop(
    "sb_anchor contains duplicated RID × assay cells (",
    nrow(duplicate_cells),
    " cells). Resolve the anchor duplication before Figure S2.",
    call. = FALSE
  )
}
sb_wide_class <- sb_class_long %>%
  select(RID, assay, class) %>%
  pivot_wider(names_from = assay, values_from = class)
missing_assay_columns <- setdiff(assay_order, names(sb_wide_class))
if (length(missing_assay_columns)) {
  stop(
    "Figure S2 is missing assay column(s): ",
    paste(missing_assay_columns, collapse = ", "),
    call. = FALSE
  )
}
sb_wide_complete <- sb_wide_class %>%
  filter(if_all(all_of(assay_order), ~ !is.na(.x)))
n_all <- nrow(sb_wide_class)
n_complete <- nrow(sb_wide_complete)
n_dropped <- n_all - n_complete
if (n_dropped > 0L) {
  message(
    "Figure S2: dropped ", n_dropped, " of ", n_all,
    " participants without all four assay classifications."
  )
}
if (n_complete == 0L) {
  stop("No participant has a complete four-assay classification panel.",
       call. = FALSE)
}
sb_wide_complete <- sb_wide_complete %>%
  unite(
    "pattern",
    all_of(assay_order),
    sep = " | ",
    remove = FALSE
  )
pattern_lookup <- sb_wide_complete %>%
  distinct(pattern, across(all_of(assay_order)))
sb_patterns <- sb_wide_complete %>%
  count(pattern, name = "n") %>%
  mutate(percent = 100 * n / n_complete) %>%
  arrange(desc(n), pattern)
n_show <- min(15L, nrow(sb_patterns))
sb_top15 <- sb_patterns %>%
  slice_head(n = n_show) %>%
  mutate(
    rank = row_number(),
    pattern_row = sprintf("%02d", rank),
    bar_label = sprintf("%d  (%.1f%%)", n, percent)
  ) %>%
  left_join(pattern_lookup, by = "pattern")
top15_coverage <- sum(sb_top15$n) / n_complete
pattern_levels <- rev(sb_top15$pattern_row)
sb_matrix <- sb_top15 %>%
  select(pattern_row, all_of(assay_order)) %>%
  pivot_longer(
    cols = all_of(assay_order),
    names_to = "assay",
    values_to = "class"
  ) %>%
  mutate(
    pattern_row = factor(pattern_row, levels = pattern_levels),
    assay = factor(assay, levels = assay_order),
    class = factor(
      class,
      levels = c("Negative", "Intermediate", "Positive")
    ),
    class_code = recode(
      as.character(class),
      "Negative" = "N",
      "Intermediate" = "I",
      "Positive" = "P"
    ),
    text_colour = if_else(class == "Positive", "white", "#263238")
  )
sb_bars <- sb_top15 %>%
  mutate(pattern_row = factor(pattern_row, levels = pattern_levels))
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
assay_labels <- setNames(assay_display(assay_order), assay_order)
FONT_FAMILY <- "sans"
COL_TEXT <- "#263238"
COL_AXIS <- "#65727A"
COL_NEGATIVE <- "#CBD8DF"
COL_INTERMEDIATE <- "#D7A044"
COL_POSITIVE <- "#315F78"
COL_BAR <- "#52636C"
COL_GRID <- "#E8ECEF"
FIG_W_MM <- 183
FIG_H_MM <- 128
row_labels <- setNames(
  paste0("#", sb_top15$rank),
  sb_top15$pattern_row
)
pattern_y_scale <- function(show_labels = TRUE) {
  scale_y_discrete(
    limits = pattern_levels,
    labels = if (show_labels) row_labels else NULL,
    drop = FALSE,
    expand = expansion(add = 0.52)
  )
}
p_matrix <- ggplot(
  sb_matrix,
  aes(x = assay, y = pattern_row, fill = class)
) +
  geom_tile(
    width = 0.90, height = 0.82,
    colour = "white", linewidth = 0.65
  ) +
  geom_text(
    aes(label = class_code, colour = text_colour),
    family = FONT_FAMILY, fontface = "bold", size = 2.35
  ) +
  scale_fill_manual(
    values = c(
      "Negative" = COL_NEGATIVE,
      "Intermediate" = COL_INTERMEDIATE,
      "Positive" = COL_POSITIVE
    ),
    breaks = c("Negative", "Intermediate", "Positive"),
    name = "Classification"
  ) +
  scale_colour_identity() +
  scale_x_discrete(
    labels = assay_labels,
    position = "top",
    expand = expansion(add = 0.10)
  ) +
  pattern_y_scale(show_labels = TRUE) +
  labs(x = NULL, y = NULL) +
  sb_theme_classic(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    text = element_text(colour = COL_TEXT),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x.top = element_text(
      size = 7.4, face = "bold", colour = COL_TEXT,
      margin = margin(b = 5)
    ),
    axis.text.y = element_text(
      size = 7, colour = COL_TEXT,
      margin = margin(r = 4)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 7.2),
    legend.text = element_text(size = 7),
    legend.key.width = grid::unit(4.5, "mm"),
    legend.key.height = grid::unit(3.2, "mm"),
    plot.margin = margin(7, 2, 5, 7)
  )
p_bars <- ggplot(sb_bars, aes(x = n, y = pattern_row)) +
  geom_col(
    width = 0.66, fill = COL_BAR,
    colour = NA
  ) +
  geom_text(
    aes(label = bar_label),
    hjust = -0.12,
    family = FONT_FAMILY, size = 2.35,
    colour = COL_TEXT
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4),
    expand = expansion(mult = c(0, 0.27)),
    position = "top"
  ) +
  pattern_y_scale(show_labels = FALSE) +
  labs(
    x = "Participants, n",
    y = NULL,
    caption = sprintf(
      "Top %d patterns represent %.1f%% of %s complete four-assay panels.",
      n_show, 100 * top15_coverage,
      format(n_complete, big.mark = ",")
    )
  ) +
  sb_theme_classic(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    text = element_text(colour = COL_TEXT),
    panel.grid.major.x = element_line(colour = COL_GRID, linewidth = 0.28),
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x.top = element_text(size = 7, colour = COL_TEXT),
    axis.title.x.top = element_text(size = 7.6, margin = margin(b = 6)),
    plot.caption = element_text(
      size = 6.7, colour = COL_AXIS,
      hjust = 1, margin = margin(t = 6)
    ),
    plot.margin = margin(7, 12, 5, 2)
  )
p <- patchwork::wrap_plots(
  p_matrix,
  p_bars,
  nrow = 1,
  widths = c(1.05, 1.25),
  guides = "collect"
) &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "horizontal",
    plot.background = element_rect(fill = "white", colour = NA)
  )
output_base <- file.path(
  study_b_dir,
  "Figure_S2_gray_zone_patterns_top15"
)
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
cat("\nFigure S2 -- top three-class patterns\n")
cat(
  "Complete four-assay panels: ", n_complete,
  "; dropped incomplete panels: ", n_dropped,
  "; top ", n_show, " coverage: ",
  sprintf("%.1f%%", 100 * top15_coverage), "\n",
  sep = ""
)
print(
  sb_top15 %>%
    select(rank, pattern, n, percent)
)
