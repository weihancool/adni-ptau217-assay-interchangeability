if (!exists("sb_theme_classic",mode = "function"))
  stop("Run the updated 01_config.R first in this same R session.",call. = FALSE)
.sb_find <- function(file) {
  if (exists("study_b_code_dir", inherits = TRUE)) {
    return(
      file.path(
        study_b_code_dir,
        file
      )
    )
  }
  env <- Sys.getenv(
    "STUDY_B_ROOT",
    unset = ""
  )
  if (
    nzchar(env) &&
    file.exists(
      file.path(
        env,
        "CODE",
        file
      )
    )
  ) {
    return(
      file.path(
        env,
        "CODE",
        file
      )
    )
  }
  d <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = FALSE
  )
  for (i in seq_len(6)) {
    p <- file.path(
      d,
      "CODE",
      file
    )
    if (file.exists(p)) {
      return(p)
    }
    if (identical(dirname(d), d)) {
      break
    }
    d <- dirname(d)
  }
  stop(
    "Cannot locate CODE/",
    file,
    call. = FALSE
  )
}
required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "purrr"
)
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]
if (length(missing_packages)) {
  stop(
    "Please install missing package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
})
if (!exists("sb_decision_long", inherits = TRUE)) {
  stop("Run 03_core_computation.R manually first.", call. = FALSE)
}
required_columns <- c(
  "RID",
  "assay",
  "amyloid_positive",
  "abs_amyloid_gap_days",
  "plasma_date",
  "log_value"
)
missing_columns <- setdiff(
  required_columns,
  names(sb_decision_long)
)
if (length(missing_columns)) {
  stop(
    "sb_decision_long is missing: ",
    paste(
      missing_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}
if (
  !exists("study_b_assays", inherits = TRUE) ||
  !length(study_b_assays)
) {
  stop(
    "study_b_assays is unavailable or empty.",
    call. = FALSE
  )
}
sb_fixed_reference <- sb_decision_long %>%
  filter(
    amyloid_positive %in% c(0, 1),
    is.finite(log_value),
    is.finite(abs_amyloid_gap_days),
    !is.na(plasma_date)
  ) %>%
  mutate(
    assay = as.character(assay)
  ) %>%
  arrange(
    RID,
    assay,
    abs_amyloid_gap_days,
    plasma_date
  ) %>%
  distinct(
    RID,
    assay,
    .keep_all = TRUE
  )
if (!nrow(sb_fixed_reference)) {
  stop(
    "No complete participant-level PET reference observations are available.",
    call. = FALSE
  )
}
missing_assays <- setdiff(
  as.character(study_b_assays),
  unique(sb_fixed_reference$assay)
)
if (length(missing_assays)) {
  stop(
    "Figure S3 has no eligible observations for assay(s): ",
    paste(
      missing_assays,
      collapse = ", "
    ),
    call. = FALSE
  )
}
sb_fixed_thresholds <- purrr::map_dfr(
  as.character(study_b_assays),
  function(a) {
    d <- sb_fixed_reference %>%
      filter(
        assay == a
      )
    status_counts <- table(
      factor(
        d$amyloid_positive,
        levels = c(0, 1)
      )
    )
    if (any(status_counts < 2L)) {
      stop(
        "Assay ",
        a,
        " requires at least two Aβ-PET-negative and two Aβ-PET-positive ",
        "participants for threshold estimation.",
        call. = FALSE
      )
    }
    mu <- mean(
      d$log_value
    )
    sig <- stats::sd(
      d$log_value
    )
    if (
      !is.finite(mu) ||
      !is.finite(sig) ||
      sig <= 0
    ) {
      stop(
        "Assay ",
        a,
        " has a non-finite or zero log-value SD.",
        call. = FALSE
      )
    }
    z <- (
      d$log_value -
        mu
    ) / sig
    th <- sb_pick_thresholds(
      d$amyloid_positive,
      z
    )
    required_rules <- c(
      "triage",
      "confirmatory",
      "youden"
    )
    if (
      !all(
        required_rules %in%
        names(th)
      ) ||
      any(
        !is.finite(
          unname(
            th[
              required_rules
            ]
          )
        )
      )
    ) {
      stop(
        "Threshold estimation failed for assay ",
        a,
        ".",
        call. = FALSE
      )
    }
    tibble(
      assay = a,
      fixed_log_mean =
        mu,
      fixed_log_sd =
        sig,
      fixed_triage_threshold_z =
        unname(
          th[["triage"]]
        ),
      fixed_confirmatory_threshold_z =
        unname(
          th[["confirmatory"]]
        ),
      fixed_youden_threshold_z =
        unname(
          th[["youden"]]
        ),
      reference_participants =
        n_distinct(
          d$RID
        ),
      amyloid_negative_n =
        sum(
          d$amyloid_positive == 0
        ),
      amyloid_positive_n =
        sum(
          d$amyloid_positive == 1
        )
    ) %>%
      mutate(
        triage_log =
          fixed_log_mean +
          fixed_triage_threshold_z *
          fixed_log_sd,
        confirmatory_log =
          fixed_log_mean +
          fixed_confirmatory_threshold_z *
          fixed_log_sd,
        youden_log =
          fixed_log_mean +
          fixed_youden_threshold_z *
          fixed_log_sd
      )
  }
)
sb_write(
  sb_fixed_thresholds %>%
    select(
      assay,
      fixed_log_mean,
      fixed_log_sd,
      fixed_triage_threshold_z,
      fixed_confirmatory_threshold_z,
      fixed_youden_threshold_z,
      reference_participants,
      amyloid_positive_n
    ),
  "Figure_S3_full_sample_fixed_thresholds.csv"
)
.audit_threshold_path <- file.path(study_b_dir, "Table_S4_fixed_threshold_audit.csv")
if (file.exists(.audit_threshold_path)) {
  .audit_thresholds <- sb_read_result(.audit_threshold_path, check.names = FALSE)
  .cols <- c("fixed_log_mean", "fixed_log_sd", "fixed_triage_threshold_z",
             "fixed_confirmatory_threshold_z", "fixed_youden_threshold_z",
             "reference_participants", "amyloid_positive_n")
  .check <- inner_join(
    sb_fixed_thresholds %>% select(assay, all_of(.cols)),
    .audit_thresholds %>% select(assay, all_of(.cols)),
    by = "assay", suffix = c("_figS3", "_audit")
  )
  .diff <- unlist(lapply(.cols, function(v) {
    abs(.check[[paste0(v, "_figS3")]] - .check[[paste0(v, "_audit")]])
  }))
  stopifnot("Table S4 and Figure S3 fixed thresholds differ" =
              nrow(.check) == length(study_b_assays) &&
              all(is.finite(.diff)) && all(.diff < 1e-12))
}
sb_threshold_lines <- sb_fixed_thresholds %>%
  select(
    assay,
    triage_log,
    youden_log,
    confirmatory_log
  ) %>%
  pivot_longer(
    cols = -assay,
    names_to =
      "rule",
    values_to =
      "log_value"
  ) %>%
  mutate(
    rule =
      sub(
        "_log$",
        "",
        rule
      ),
    rule = factor(
      rule,
      levels = c(
        "triage",
        "youden",
        "confirmatory"
      ),
      labels = c(
        "Triage",
        "Youden",
        "Confirmatory"
      )
    ),
    threshold_tag = dplyr::recode(
      as.character(rule),
      "Triage" = "T",
      "Youden" = "Y",
      "Confirmatory" = "C"
    ),
    threshold_vjust = dplyr::recode(
      as.character(rule),
      "Triage" = -1.2,
      "Youden" = -2.8,
      "Confirmatory" = -4.4
    )
  )
unexpected_order <- sb_fixed_thresholds %>%
  filter(
    !(
      triage_log <=
        youden_log &
        youden_log <=
        confirmatory_log
    )
  )
if (nrow(unexpected_order)) {
  warning(
    "Threshold order is not triage ≤ Youden ≤ confirmatory for: ",
    paste(
      unexpected_order$assay,
      collapse = ", "
    ),
    ". Lines are plotted at their estimated values.",
    call. = FALSE
  )
}
assay_display <- function(x) {
  dplyr::recode(
    as.character(x),
    "AlzPath" =
      "ALZpath",
    "ALZpath" =
      "ALZpath",
    "C2N_ratio" =
      "C2N",
    "C2N" =
      "C2N",
    .default =
      as.character(x)
  )
}
sb_facet_labels <- sb_fixed_thresholds %>%
  transmute(
    assay,
    assay_label =
      assay_display(
        assay
      )
  )
counts_are_common <-
  n_distinct(
    sb_fixed_thresholds$amyloid_negative_n
  ) == 1L &&
  n_distinct(
    sb_fixed_thresholds$amyloid_positive_n
  ) == 1L
if (counts_are_common) {
  n_negative <- unique(
    sb_fixed_thresholds$amyloid_negative_n
  )
  n_positive <- unique(
    sb_fixed_thresholds$amyloid_positive_n
  )
  sb_facet_labels <- sb_facet_labels %>%
    mutate(
      facet_label =
        assay_label
    )
  pet_legend_labels <- c(
    "Aβ-PET negative" =
      sprintf(
        "Aβ-PET−  (n = %d)",
        n_negative
      ),
    "Aβ-PET positive" =
      sprintf(
        "Aβ-PET+  (n = %d)",
        n_positive
      )
  )
} else {
  sb_facet_labels <- sb_facet_labels %>%
    left_join(
      select(
        sb_fixed_thresholds,
        assay,
        amyloid_negative_n,
        amyloid_positive_n
      ),
      by = "assay"
    ) %>%
    mutate(
      facet_label =
        sprintf(
          "%s\nAβ− n = %d   Aβ+ n = %d",
          assay_label,
          amyloid_negative_n,
          amyloid_positive_n
        )
    )
  pet_legend_labels <- c(
    "Aβ-PET negative" =
      "Aβ-PET−",
    "Aβ-PET positive" =
      "Aβ-PET+"
  )
}
facet_levels <- sb_facet_labels$facet_label[
  match(
    as.character(
      study_b_assays
    ),
    sb_facet_labels$assay
  )
]
sb_plot_data <- sb_fixed_reference %>%
  mutate(
    pet_status =
      factor(
        if_else(
          amyloid_positive == 1,
          "Aβ-PET positive",
          "Aβ-PET negative"
        ),
        levels = c(
          "Aβ-PET negative",
          "Aβ-PET positive"
        )
      )
  ) %>%
  left_join(
    select(
      sb_facet_labels,
      assay,
      facet_label
    ),
    by = "assay"
  ) %>%
  mutate(
    facet_label =
      factor(
        facet_label,
        levels = facet_levels
      )
  )
sb_threshold_lines <- sb_threshold_lines %>%
  left_join(
    select(
      sb_facet_labels,
      assay,
      facet_label
    ),
    by = "assay"
  ) %>%
  mutate(
    facet_label =
      factor(
        facet_label,
        levels = facet_levels
      )
  )
FONT_FAMILY <- "sans"
COL_TEXT <- "#263238"
COL_AXIS <- "#65727A"
COL_NEG <- "#5E879E"
COL_POS <- "#D49A45"
COL_THRESHOLD <- "#353B40"
FIG_W_MM <- 183
FIG_H_MM <- 114
p <- ggplot(
  sb_plot_data,
  aes(
    x = log_value,
    y = pet_status
  )
) +
geom_violin(
  aes(
    fill = pet_status,
    colour = pet_status
  ),
  width = 0.80,
  scale = "width",
  trim = TRUE,
  alpha = 0.25,
  linewidth = 0.45,
  show.legend = TRUE,
  key_glyph = "rect"
) +
geom_point(
  aes(
    colour = pet_status
  ),
  position =
    position_jitter(
      width = 0,
      height = 0.065,
      seed = 20260910L
    ),
  shape = 16,
  size = 0.50,
  alpha = 0.18,
  show.legend = FALSE
) +
geom_boxplot(
  width = 0.13,
  outlier.shape = NA,
  fill = "white",
  colour = COL_TEXT,
  linewidth = 0.38,
  alpha = 0.88,
  show.legend = FALSE
) +
geom_vline(
  data =
    sb_threshold_lines,
  aes(
    xintercept = log_value,
    linetype = rule
  ),
  colour =
    COL_THRESHOLD,
  linewidth =
    0.58,
  show.legend =
    TRUE,
  key_glyph =
    "path",
  inherit.aes =
    FALSE
) +
  geom_text(
    data = sb_threshold_lines,
    aes(
      x = log_value,
      y = "Aβ-PET positive",
      label = threshold_tag,
      vjust = threshold_vjust
    ),
    colour = COL_THRESHOLD,
    size = 2.15,
    fontface = "bold",
    family = FONT_FAMILY,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  scale_y_discrete(
    expand = expansion(add = c(0.28, 0.78))
  ) +
facet_wrap(
  ~ facet_label,
  ncol = 2,
  scales = "free_x"
) +
scale_fill_manual(
  values = c(
    "Aβ-PET negative" =
      COL_NEG,
    "Aβ-PET positive" =
      COL_POS
  ),
  breaks = c(
    "Aβ-PET negative",
    "Aβ-PET positive"
  ),
  labels =
    pet_legend_labels,
  name =
    "Amyloid PET status"
) +
  scale_colour_manual(
    values = c(
      "Aβ-PET negative" =
        COL_NEG,
      "Aβ-PET positive" =
        COL_POS
    ),
    guide =
      "none"
  ) +
scale_linetype_manual(
  values = c(
    "Triage" =
      "dashed",
    "Youden" =
      "solid",
    "Confirmatory" =
      "dotdash"
  ),
  breaks = c(
    "Triage",
    "Youden",
    "Confirmatory"
  ),
  name =
    "Fixed threshold"
) +
labs(
  x =
    "Log-transformed assay value (assay-specific scale)",
  y =
    NULL
) +
guides(
  fill = guide_legend(
    order = 1,
    nrow = 1,
    byrow = TRUE,
    title.position = "left",
    title.hjust = 0,
    label.position = "right",
    keywidth =
      grid::unit(
        6.0,
        "mm"
      ),
    keyheight =
      grid::unit(
        3.0,
        "mm"
      ),
    override.aes = list(
      alpha = 0.48,
      linewidth = 0.45
    )
  ),
  linetype = guide_legend(
    order = 2,
    nrow = 1,
    byrow = TRUE,
    title.position = "left",
    title.hjust = 0,
    label.position = "right",
    keywidth =
      grid::unit(
        9.0,
        "mm"
      ),
    keyheight =
      grid::unit(
        2.8,
        "mm"
      ),
    override.aes = list(
      colour =
        COL_THRESHOLD,
      linewidth =
        0.85
    )
  )
) +
sb_theme_classic(
  base_size = 8,
  base_family = FONT_FAMILY
) +
  theme(
    text =
      element_text(
        colour = COL_TEXT
      ),
    axis.title.x =
      element_text(
        size = 8,
        margin = margin(
          t = 7
        )
      ),
    axis.text.x =
      element_text(
        size = 7.2,
        colour = COL_TEXT
      ),
    axis.text.y =
      element_text(
        size = 7.4,
        colour = COL_TEXT,
        margin = margin(
          r = 5
        )
      ),
    axis.line =
      element_line(
        colour = COL_AXIS,
        linewidth = 0.32
      ),
    axis.ticks =
      element_line(
        colour = COL_AXIS,
        linewidth = 0.32
      ),
    axis.ticks.length =
      grid::unit(
        1.2,
        "mm"
      ),
    strip.background =
      element_blank(),
    strip.text =
      element_text(
        size = 8.2,
        face = "bold",
        lineheight = 0.98,
        margin = margin(
          b = 4
        )
      ),
    panel.spacing.x =
      grid::unit(
        9,
        "mm"
      ),
    panel.spacing.y =
      grid::unit(
        6,
        "mm"
      ),
    legend.position =
      "bottom",
    legend.box =
      "horizontal",
    legend.box.just =
      "center",
    legend.justification =
      "center",
    legend.title =
      element_text(
        size = 7.4,
        face = "bold",
        colour = COL_TEXT,
        hjust = 0,
        margin = margin(
          r = 4.0,
          b = 0,
          unit = "pt"
        )
      ),
    legend.text =
      element_text(
        size = 7.2,
        colour = COL_TEXT,
        margin = margin(
          l = 1,
          r = 6,
          unit = "pt"
        )
      ),
    legend.key =
      element_blank(),
    legend.spacing.x =
      grid::unit(
        1.8,
        "mm"
      ),
    legend.spacing.y =
      grid::unit(
        0.2,
        "mm"
      ),
    legend.margin =
      margin(
        t = 2,
        b = 0,
        unit = "pt"
      ),
    legend.box.spacing =
      grid::unit(
        1,
        "pt"
      ),
    plot.background =
      element_rect(
        fill = "white",
        colour = NA
      ),
    plot.margin =
      margin(
        8,
        9,
        5,
        8
      )
  )
output_base <- file.path(
  study_b_dir,
  "Figure_S3_threshold_distributions"
)
pdf_device <- if (
  capabilities("cairo")
) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}
sb_ggsave(
  paste0(
    output_base,
    ".pdf"
  ),
  p,
  width =
    FIG_W_MM,
  height =
    FIG_H_MM,
  units =
    "mm",
  device =
    pdf_device,
  bg =
    "white"
)
sb_ggsave(
  paste0(
    output_base,
    ".tiff"
  ),
  p,
  width =
    FIG_W_MM,
  height =
    FIG_H_MM,
  units =
    "mm",
  dpi =
    600,
  compression =
    "lzw",
  bg =
    "white"
)
sb_ggsave(
  paste0(
    output_base,
    ".png"
  ),
  p,
  width =
    FIG_W_MM,
  height =
    FIG_H_MM,
  units =
    "mm",
  dpi =
    600,
  bg =
    "white"
)
cat(
  "\nFigure S3 -- full-sample fixed thresholds\n"
)
print(
  sb_fixed_thresholds %>%
    select(
      assay,
      fixed_log_mean,
      fixed_log_sd,
      fixed_triage_threshold_z,
      fixed_confirmatory_threshold_z,
      fixed_youden_threshold_z,
      reference_participants,
      amyloid_positive_n
    )
)
.ref_counts <- sb_fixed_reference %>%
  distinct(RID, amyloid_positive)
if (anyDuplicated(.ref_counts$RID))
  stop("Figure S3 reference PET status differs across assays for one RID.", call. = FALSE)
writeLines(paste0(
  "Figure S3. Assay distributions and full-sample fixed thresholds. ",
  "One eligible PET-matched plasma collection was selected per participant ",
  "by the smallest absolute plasma-PET interval; ties were resolved by the earlier plasma collection. ",
  "The reference set includes ", nrow(.ref_counts), " participants (",
  sum(.ref_counts$amyloid_positive == 0L), " PET negative and ",
  sum(.ref_counts$amyloid_positive == 1L), " PET positive) and differs from the primary-anchor PET evaluation set. ",
  "These in-sample fixed thresholds support the calibration sensitivity analyses ",
  "and do not replace the primary cross-fitted thresholds."
), file.path(study_b_dir, "Figure_S3_caption.txt"))
