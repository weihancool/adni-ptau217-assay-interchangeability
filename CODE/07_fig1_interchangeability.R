if (!exists("sb_theme_classic",mode = "function"))
  stop("Run the updated 01_config.R first in this same R session.",call. = FALSE)
pkgs <- c("ggplot2", "dplyr", "patchwork", "tidyr", "purrr", "tibble")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Please install missing package(s): ",
       paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(purrr)
library(tibble)
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
if (exists("study_b_dir", inherits = TRUE)) {
  RESULTS_DIR <- study_b_dir
} else {
  project_root <- Sys.getenv("STUDY_B_ROOT", unset = "")
  if (!nzchar(project_root)) {
    stop("Please set STUDY_B_ROOT, for example:\n",
         'Source CODE/01_config.R first or set STUDY_B_ROOT to this project root.', call. = FALSE)
  }
  RESULTS_DIR <- file.path(project_root, "RESULTS")
}
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
file_A     <- file.path(RESULTS_DIR, "Figure_1A_AARR.csv")
file_C     <- file.path(RESULTS_DIR, "Figure_1C_three_class_composition.csv")
file_panel <- file.path(RESULTS_DIR, "Figure_1_participant_panel.csv")
if (!exists("sb_assert_code_current", mode = "function"))
  stop("Run the updated 01_config.R first.", call. = FALSE)
sb_assert_code_current()
if (!exists("sb_core_contract_used", inherits = TRUE) ||
    !identical(sb_core_contract_used, sb_make_core_contract()))
  stop("Core objects are stale. Run 02_preprocess.R and 03_core_computation.R ",
       "after the current 01_config.R in a clean session.", call. = FALSE)
if (!sb_have("sb_anchor", "sb_events", "sb_crossfit"))
  stop("Run 02_preprocess.R and 03_core_computation.R manually first.", call. = FALSE)
sb_aarr_point <- vapply(sb_rules, function(col)
  sb_discordant_rate(sb_rule_matrix(sb_anchor, col)), numeric(1))
sb_aarr_conditional <- lapply(sb_rules, function(col)
  sb_boot_summary(sb_anchor, sb_aarr_stat(col), B = study_b_boot_n,
                  seed = sb_seed(paste0("aarr_conditional_", col))))
.sb_ev_index <- split(seq_len(nrow(sb_events)), as.character(sb_events$RID))
.sb_ev_ids   <- names(.sb_ev_index)
.sb_pipeline_one <- function(b) {
  set.seed(study_b_seed_aarr + b)
  draw <- sample(.sb_ev_ids, length(.sb_ev_ids), replace = TRUE)
  d <- bind_rows(lapply(seq_along(draw), function(j) {
    z <- sb_events[.sb_ev_index[[draw[j]]], , drop = FALSE]
    z$original_RID <- z$RID
    z$RID <- j
    z
  }))
  d$sample_id <- seq_len(nrow(d))
  d <- sb_assign_boot_folds(d, seed = study_b_seed_aarr + 100000L + b)
  sc <- sb_crossfit(d)
  anch <- filter(sc, primary_anchor)
  vapply(sb_rules, function(col) sb_discordant_rate(sb_rule_matrix(anch, col)),
         numeric(1))
}
sb_aarr_pipeline_draws <- matrix(NA_real_, nrow = study_b_boot_n_pipeline,
                                 ncol = length(sb_rules),
                                 dimnames = list(NULL, names(sb_rules)))
sb_aarr_pipeline_errors <- rep(NA_character_, study_b_boot_n_pipeline)
if (study_b_boot_n_pipeline > 0) {
  message("Figure 1: full-pipeline bootstrap, ", study_b_boot_n_pipeline, " replicates.")
  .every <- max(1L, floor(study_b_boot_n_pipeline / 20L))
  for (b in seq_len(study_b_boot_n_pipeline)) {
    sb_aarr_pipeline_draws[b, ] <- tryCatch(.sb_pipeline_one(b), error = function(e) {
      sb_aarr_pipeline_errors[b] <<- conditionMessage(e)
      rep(NA_real_, length(sb_rules))
    })
    if (b %% .every == 0L) message("  pipeline bootstrap ", b, "/", study_b_boot_n_pipeline)
  }
}
sb_write(data.frame(draw = seq_len(nrow(sb_aarr_pipeline_draws)), sb_aarr_pipeline_draws,
                    error = sb_aarr_pipeline_errors),
         "Figure_1A_AARR_pipeline_bootstrap_draws.csv")
sb_write(data.frame(draw = which(!is.na(sb_aarr_pipeline_errors)),
                    error = sb_aarr_pipeline_errors[!is.na(sb_aarr_pipeline_errors)]),
         "Figure_1A_AARR_pipeline_bootstrap_failures.csv")
message("Figure 1 pipeline bootstrap: ",
        sum(complete.cases(sb_aarr_pipeline_draws)), "/", study_b_boot_n_pipeline,
        " complete replicates retained.")
if (any(!is.na(sb_aarr_pipeline_errors))) {
  cat("\nFigure 1 pipeline-bootstrap failures by reason:\n")
  print(table(sb_aarr_pipeline_errors, useNA = "no"))
}
sb_overall <- tibble(
  rule = names(sb_rules),
  n_participants = vapply(sb_rules, function(col)
    nrow(sb_complete_panel(sb_rule_matrix(sb_anchor, col))), numeric(1)),
  AARR = unname(sb_aarr_point),
  AARR_low  = vapply(names(sb_rules), function(r) sb_ci(sb_aarr_pipeline_draws[, r])[["low"]],  numeric(1)),
  AARR_high = vapply(names(sb_rules), function(r) sb_ci(sb_aarr_pipeline_draws[, r])[["high"]], numeric(1)),
  pipeline_bootstrap_successful = colSums(is.finite(sb_aarr_pipeline_draws)),
  pipeline_bootstrap_attempted = study_b_boot_n_pipeline,
  conditional_low  = vapply(sb_aarr_conditional, function(z) unname(z[["low"]]),  numeric(1)),
  conditional_high = vapply(sb_aarr_conditional, function(z) unname(z[["high"]]), numeric(1)),
  CI_scope = "Primary CI: participant bootstrap re-running folds, standardisation, thresholds, calibration and classification"
)
if (any(sb_overall$pipeline_bootstrap_successful < 0.90 * study_b_boot_n_pipeline))
  warning("Figure 1 full-pipeline bootstrap success below 90%; inspect the failure audit.", call. = FALSE)
sb_write(sb_overall, "Figure_1A_AARR.csv")
sb_pairwise <- imap_dfr(sb_rules, function(col, rule) {
  map_dfr(sb_pairs, function(p) {
    d <- filter(sb_anchor, assay %in% p)
    w <- d %>% select(RID, assay, value = all_of(col)) %>%
      pivot_wider(names_from = assay, values_from = value)
    z <- sb_agreement(w[[p[1]]], w[[p[2]]])
    ag <- sb_boot_summary(d, sb_pair_stat(col, p, "agreement"),
                          B = study_b_boot_n, seed = study_b_seed_pair)
    kp <- sb_boot_summary(d, sb_pair_stat(col, p, "kappa"),
                          B = study_b_boot_n, seed = study_b_seed_pair + 1L)
    ac <- sb_boot_summary(d, sb_pair_stat(col, p, "gwet_ac1"),
                          B = study_b_boot_n, seed = study_b_seed_pair + 2L)
    tibble(rule = rule, assay_1 = p[1], assay_2 = p[2], n = nrow(w),
           agreement = unname(z["agreement"]),
           agreement_low = unname(ag["low"]), agreement_high = unname(ag["high"]),
           kappa = unname(z["kappa"]),
           kappa_low = unname(kp["low"]), kappa_high = unname(kp["high"]),
           gwet_ac1 = unname(z["gwet_ac1"]),
           gwet_ac1_low = unname(ac["low"]), gwet_ac1_high = unname(ac["high"]))
  })
})
sb_write(sb_pairwise, "Figure_1B_pairwise_agreement_all_rules.csv")
.sb_three_class_shares <- function(d) {
  id <- if (".boot_RID" %in% names(d)) ".boot_RID" else "RID"
  w <- d %>% select(id = all_of(id), assay, three_class) %>%
    mutate(three_class = as.character(three_class)) %>%
    pivot_wider(names_from = assay, values_from = three_class)
  m <- sb_complete_panel(as.matrix(w[, study_b_assays]))
  comp <- apply(m, 1, function(z) length(unique(z)) == 1L)
  ext  <- apply(m, 1, function(z) all(c("negative", "positive") %in% z))
  c(complete = mean(comp), adjacent = mean(!comp & !ext), extreme = mean(ext))
}
sb_gray_point <- .sb_three_class_shares(sb_anchor)
sb_gray_draws <- sb_cluster_boot_multi(sb_anchor, .sb_three_class_shares,
                                       B = study_b_boot_n, seed = study_b_seed_gray)
sb_fig1c <- tibble(
  category = c("Complete agreement", "Gray-zone / adjacent only", "Extreme negative-positive"),
  n = as.integer(round(sb_gray_point * n_distinct(sb_anchor$RID))),
  percent = 100 * unname(sb_gray_point),
  low  = 100 * vapply(colnames(sb_gray_draws), function(k) sb_ci(sb_gray_draws[, k])[["low"]],  numeric(1)),
  high = 100 * vapply(colnames(sb_gray_draws), function(k) sb_ci(sb_gray_draws[, k])[["high"]], numeric(1)),
  bootstrap_successful = colSums(is.finite(sb_gray_draws)),
  bootstrap_attempted = study_b_boot_n,
  CI_scope = "Participant bootstrap conditional on frozen out-of-fold predictions"
)
sb_write(sb_fig1c, "Figure_1C_three_class_composition.csv")
cat("\nFigure 1 -- cross-assay disagreement\n"); print(sb_overall)
cat("\nFigure 1 -- three-level composition\n"); print(sb_fig1c)
sb_fig1_auc <- purrr::map_dfr(study_b_assays, function(a) {
  g <- sb_anchor %>% filter(assay == a, !is.na(amyloid_positive))
  if (!nrow(g)) stop("No PET reference for assay ", a, call. = FALSE)
  ci <- sb_boot_summary(
    g, function(d) sb_auc_rank(d$amyloid_positive, d$calibrated_probability),
    B = study_b_boot_n, seed = sb_seed(paste0("auc_", a)))
  tibble(assay = a, n = nrow(g),
         crossfit_AUC = sb_auc_rank(g$amyloid_positive, g$calibrated_probability),
         AUC_low = unname(ci["low"]), AUC_high = unname(ci["high"]),
         AUC_bootstrap_successful = unname(ci["bootstrap_successful"]),
         bootstrap_attempted = study_b_boot_n,
         CI_scope = "Participant bootstrap conditional on frozen out-of-fold predictions")
})
if (any(sb_fig1_auc$AUC_bootstrap_successful < 0.90 * study_b_boot_n))
  warning("Figure 1 AUC bootstrap success below 90%; report the retained counts.", call. = FALSE)
sb_write(sb_fig1_auc, "Figure_1_AUC.csv")
assay_label <- function(x) {
  dplyr::recode(trimws(as.character(x)),
                "AlzPath" = "ALZpath", "ALZpath" = "ALZpath",
                "C2N_ratio" = "C2N", "C2N" = "C2N",
                .default = trimws(as.character(x)))
}
.pattern_of <- function(z) {
  if (length(unique(z)) == 1L) return("Complete agreement")
  if (all(c("negative", "positive") %in% z)) return("Extreme conflict")
  "Intermediate / adjacent"
}
pattern_levels <- c("Complete agreement", "Intermediate / adjacent",
                    "Extreme conflict")
assay_levels <- assay_label(study_b_assays)
class_wide <- sb_anchor %>%
  select(RID, assay, three_class) %>%
  mutate(three_class = as.character(three_class)) %>%
  tidyr::pivot_wider(names_from = assay, values_from = three_class)
missing_assay_cols <- setdiff(study_b_assays, names(class_wide))
if (length(missing_assay_cols)) {
  stop("Figure 1: sb_anchor is missing assay column(s): ",
       paste(missing_assay_cols, collapse = ", "), call. = FALSE)
}
class_mat <- as.matrix(class_wide[, study_b_assays])
if (any(rowSums(is.na(class_mat)) > 0)) {
  stop("Figure 1: sb_anchor contains an incomplete four-assay panel; ",
       "panels A-D must be drawn on complete panels only.", call. = FALSE)
}
participant_meta <- tibble(
  RID = class_wide$RID,
  pattern = apply(class_mat, 1, .pattern_of)
) %>%
  left_join(
    sb_anchor %>%
      group_by(RID) %>%
      summarise(mean_margin = mean(margin_youden),
                amyloid_positive = sb_first_nonmissing(amyloid_positive),
                .groups = "drop"),
    by = "RID"
  ) %>%
  arrange(desc(mean_margin)) %>%
  mutate(
    rank = row_number(),
    pet_key = dplyr::case_when(
      is.na(amyloid_positive)   ~ "No PET reference",
      amyloid_positive == 1     ~ "Amyloid PET positive",
      TRUE                      ~ "Amyloid PET negative"),
    pattern = factor(pattern, levels = pattern_levels)
  )
participant_panel <- sb_anchor %>%
  transmute(RID, assay,
            assay_label = factor(assay_label(assay), levels = assay_levels),
            three_class = as.character(three_class),
            margin_youden) %>%
  left_join(select(participant_meta, RID, rank, pattern, pet_key, mean_margin),
            by = "RID")
sb_write(
  participant_panel %>%
    mutate(assay_label = as.character(assay_label),
           pattern = as.character(pattern)) %>%
    arrange(rank, assay_label),
  "Figure_1_participant_panel.csv"
)
participant_meta <- participant_meta %>%
  select(RID, rank, pattern, pet_key, mean_margin) %>%
  arrange(rank)
local({
  results <- RESULTS_DIR
  inputs <- c("Figure_1A_AARR.csv", "Figure_1C_three_class_composition.csv",
              "Figure_1_participant_panel.csv", "Figure_1_AUC.csv")
  paths <- file.path(results, inputs)
  read_input <- function(path, columns) {
    if (!file.exists(path)) stop("Missing newly computed Figure 1 input: ", path,
                                 call. = FALSE)
    z <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(z) || !all(columns %in% names(z)))
      stop("Empty or incompatible newly computed Figure 1 input: ", path,
           call. = FALSE)
    z
  }
  rules <- read_input(paths[1], c("rule", "AARR", "AARR_low", "AARR_high"))
  composition <- read_input(paths[2], c("category", "percent", "low", "high"))
  panel <- read_input(paths[3], c("RID", "assay_label", "margin_youden", "rank", "pattern"))
  performance <- read_input(paths[4], c("assay", "crossfit_AUC", "AUC_low", "AUC_high"))
  panel$assay_label <- trimws(as.character(panel$assay_label))
  rules$rule <- trimws(as.character(rules$rule))
  rule_order <- c("triage", "confirmatory", "youden")
  if (nrow(rules) != 3L || !setequal(rules$rule, rule_order))
    stop("AARR input must contain one row for each of triage, confirmatory and youden.", call. = FALSE)
  numeric_rules <- unlist(rules[c("AARR", "AARR_low", "AARR_high")])
  if (!is.numeric(numeric_rules) || any(!is.finite(numeric_rules)) ||
      any(numeric_rules < 0 | numeric_rules > 1) || any(rules$AARR_low > rules$AARR_high))
    stop("Computed AARR estimates/intervals must be finite proportions in [0,1].", call. = FALSE)
  if (!is.numeric(panel$margin_youden) || any(!is.finite(panel$margin_youden)) ||
      anyNA(panel$RID) || anyNA(panel$assay_label) ||
      anyDuplicated(panel[c("RID", "assay_label")]) ||
      any(table(panel$RID) != 4L) || length(unique(panel$assay_label)) != 4L)
    stop("Participant input must contain exactly four finite assay margins per participant.", call. = FALSE)
  spans <- panel %>% group_by(RID) %>%
    summarise(lo = min(margin_youden), hi = max(margin_youden),
              rank = first(rank), pattern = first(pattern),
              n_rank = n_distinct(rank), n_pattern = n_distinct(pattern), .groups = "drop") %>%
    arrange(rank) %>% mutate(discordant = lo < 0 & hi >= 0)
  n <- nrow(spans)
  if (any(spans$n_rank != 1L | spans$n_pattern != 1L) ||
      anyNA(spans$rank) || !identical(as.numeric(spans$rank), as.numeric(seq_len(n))))
    stop("Each participant must have one unique, consecutive plotting rank and pattern.", call. = FALSE)
  k <- sum(spans$discordant)
  youden <- rules[rules$rule == "youden", , drop = FALSE]
  if (abs(k / n - youden$AARR) > 1e-6)
    stop("Computed participant margins and Youden AARR disagree.", call. = FALSE)
  if ("n_participants" %in% names(rules) && any(rules$n_participants != n))
    stop("AARR and participant inputs have different cohort sizes.", call. = FALSE)
  category_map <- c("Complete agreement" = "All four agree",
                    "Gray-zone / adjacent only" = "Adjacent-only disagreement",
                    "Intermediate / adjacent" = "Adjacent-only disagreement",
                    "Extreme negative-positive" = "Negative-positive conflict",
                    "Extreme conflict" = "Negative-positive conflict")
  composition$label <- unname(category_map[trimws(as.character(composition$category))])
  category_order <- c("All four agree", "Adjacent-only disagreement", "Negative-positive conflict")
  if (nrow(composition) != 3L || anyNA(composition$label) || anyDuplicated(composition$label) ||
      !is.numeric(composition$percent) || any(!is.finite(composition$percent)) ||
      any(composition$percent < 0 | composition$percent > 100) ||
      abs(sum(composition$percent) - 100) > 1e-6)
    stop("Three-level input must contain three mutually exclusive categories summing to 100%.", call. = FALSE)
  observed <- c("All four agree" = sum(spans$pattern == "Complete agreement"),
                "Adjacent-only disagreement" = sum(spans$pattern == "Intermediate / adjacent"),
                "Negative-positive conflict" = sum(spans$pattern == "Extreme conflict"))
  if (sum(observed) != n ||
      any(abs(composition$percent - 100 * observed[composition$label] / n) > 1e-6))
    stop("Frozen three-level composition and participant input disagree.", call. = FALSE)
  ink <- "#253444"; blue <- "#254F70"; red <- "#B44145"
  quiet <- "#BBC5CC"; adjacent <- "#8597A6"; grid_col <- "#E7EBEE"
  font <- "sans"
  base_theme <- sb_theme_classic() + theme(legend.position = "none")
  lim <- range(panel$margin_youden)
  lim <- c(floor(min(lim[1], 0) * 2) / 2, ceiling(max(lim[2], 0) * 2) / 2)
  if (diff(lim) == 0) lim <- lim + c(-0.5, 0.5)
  pB <- ggplot() +
    annotate("rect",xmin = lim[1],xmax = 0,ymin = 0,ymax = lim[2],
             fill = "#F7EBEB",colour = NA) +
    geom_abline(slope = 1,intercept = 0,colour = grid_col,
                linewidth = 0.35,linetype = "dashed") +
    geom_vline(xintercept = 0, colour = ink, linewidth = 0.45) +
    geom_hline(yintercept = 0, colour = ink, linewidth = 0.45) +
    geom_point(data = filter(spans, !discordant),aes(x = lo,y = hi),
               colour = adjacent,size = 0.85,alpha = 0.50) +
    geom_point(data = filter(spans, discordant),aes(x = lo,y = hi),
               colour = red,size = 1.1,alpha = 0.82) +
    scale_x_continuous(breaks = pretty(lim, n = 5), expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = pretty(lim, n = 5), expand = expansion(mult = 0)) +
    coord_equal(xlim = lim,ylim = lim,expand = FALSE) +
    labs(title = "B",
         x = "Smallest Youden margin (SD units)",
         y = "Largest Youden margin (SD units)") +
    base_theme + theme(axis.title.x = element_text(margin = margin(t = 6)))
  normalise_assay <- function(x) dplyr::recode(trimws(as.character(x)),
                                               "AlzPath" = "ALZpath", "ALZpath" = "ALZpath", "C2N_ratio" = "C2N", .default = trimws(as.character(x)))
  performance$assay <- normalise_assay(performance$assay)
  auc <- performance %>% select(assay,crossfit_AUC,AUC_low,AUC_high) %>% distinct()
  assay_order <- c("ALZpath","C2N","Fujirebio","Janssen")
  if (nrow(auc) != 4L || !setequal(auc$assay,assay_order))
    stop("AUC input requires one consistent estimate and CI per assay across decision rules.",call. = FALSE)
  z <- unlist(auc[c("crossfit_AUC","AUC_low","AUC_high")])
  if (!is.numeric(z) || any(!is.finite(z)) || any(z < 0 | z > 1) || any(auc$AUC_low > auc$AUC_high))
    stop("Invalid frozen AUC estimates or intervals.",call. = FALSE)
  auc <- auc[match(assay_order,auc$assay),]
  auc$y <- 4:1
  auc$annotation <- sprintf("%.3f [%.3f, %.3f]",auc$crossfit_AUC,auc$AUC_low,auc$AUC_high)
  auc_floor <- min(0.75, floor(min(auc$AUC_low) * 20) / 20)
  pet_n <- if ("n" %in% names(performance)) unique(performance$n) else numeric()
  if (length(pet_n) && (any(!is.finite(pet_n)) || any(pet_n < 1 | pet_n > n)))
    stop("PET performance and Figure 1 cohort sizes are incompatible.",call. = FALSE)
  pA <- ggplot(auc) +
    geom_vline(xintercept = seq(auc_floor,1,0.05),colour = grid_col,linewidth = 0.3) +
    geom_segment(aes(x = AUC_low,xend = AUC_high,y = y,yend = y),colour = blue,linewidth = 0.8) +
    geom_point(aes(x = crossfit_AUC,y = y),colour = blue,size = 2.3) +
    geom_text(aes(x = 1.018,y = y,label = annotation),
              hjust = 0,size = 2.65,family = font,colour = ink) +
    scale_x_continuous(limits = c(auc_floor,1.20),breaks = seq(auc_floor,1,0.05),
                       expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = auc$y,labels = auc$assay,limits = c(0.55,4.45),
                       expand = expansion(mult = 0)) +
    labs(title = "A",x = "AUC (95% CI)",y = NULL) + base_theme +
    theme(axis.line.y = element_blank(),axis.ticks.y = element_blank())
  binary <- rules[match(rule_order,rules$rule),]
  row_names <- c("Triage","Confirmatory","Youden")
  binary$y <- 3:1
  binary_caption <- paste(sprintf("%s %.1f%% [%.1f-%.1f]",row_names,100*binary$AARR,
                                  100*binary$AARR_low,100*binary$AARR_high),collapse = "; ")
  pC <- ggplot(binary,aes(y = y)) +
    geom_vline(xintercept = seq(0,100,25),colour = grid_col,linewidth = 0.3) +
    geom_segment(aes(x = 100*AARR_low,xend = 100*AARR_high,yend = y),
                 colour = blue,linewidth = 0.8) +
    geom_point(aes(x = 100*AARR),colour = blue,size = 2.3) +
    geom_text(aes(x = 102,y = y,
                  label = sprintf("%.1f [%.1f, %.1f]",100*AARR,100*AARR_low,100*AARR_high)),
              hjust = 0,size = 2.5,family = font,colour = ink) +
    scale_x_continuous(limits = c(0,151),breaks = seq(0,100,25),expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = 3:1,labels = row_names,limits = c(0.5,3.5),expand = expansion(mult = 0)) +
    labs(title = "C",x = "Four-assay disagreement (%)",y = NULL) + base_theme +
    theme(axis.line.y = element_blank(),axis.ticks.y = element_blank())
  composition <- composition[match(category_order,composition$label),]
  if (any(!is.finite(unlist(composition[c("low","high")]))) ||
      any(composition$low < 0 | composition$high > 100 | composition$low > composition$high))
    stop("Invalid three-level composition intervals.",call. = FALSE)
  composition$y <- 3:1
  conflict <- composition[composition$label == "Negative-positive conflict",]
  conflict_n <- observed[["Negative-positive conflict"]]
  pD <- ggplot(composition,aes(y = y,colour = label)) +
    geom_vline(xintercept = seq(0,100,25),colour = grid_col,linewidth = 0.3) +
    geom_segment(aes(x = low,xend = high,yend = y),linewidth = 0.8) +
    geom_point(aes(x = percent),size = 2.3) +
    geom_text(aes(x = 102,label = sprintf("%.1f%%",percent)),hjust = 0,
              size = 2.65,family = font) +
    scale_colour_manual(values = c("All four agree" = adjacent,
                                   "Adjacent-only disagreement" = blue,
                                   "Negative-positive conflict" = red),guide = "none") +
    scale_x_continuous(limits = c(0,122),breaks = seq(0,100,25),expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = 3:1,
                       labels = c("All four agree","Adjacent-only\ndisagreement","Negative-positive\nconflict"),
                       limits = c(0.5,3.5),expand = expansion(mult = 0)) +
    labs(title = "D",x = "Three-level classification patterns (%)",y = NULL) + base_theme +
    theme(axis.line.y = element_blank(),axis.ticks.y = element_blank())
  fig <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1.15,0.85)) +
    plot_annotation(theme = theme(plot.background = element_rect(fill = "white",colour = NA),
                                  plot.margin = margin(2,2,2,2)))
  out <- file.path(results,"Figure_1_interchangeability")
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  sb_ggsave(paste0(out,".pdf"),fig,width = 190,height = 150,units = "mm",device = pdf_device,bg = "white")
  sb_ggsave(paste0(out,".png"),fig,width = 190,height = 150,units = "mm",dpi = 600,bg = "white")
  sb_ggsave(paste0(out,".tiff"),fig,width = 190,height = 150,units = "mm",dpi = 600,compression = "lzw",bg = "white")
  sb_write(bind_rows(auc %>% mutate(panel = "A"),spans %>% mutate(panel = "B"),
                     binary %>% mutate(panel = "C"),composition %>% mutate(panel = "D")),
           "Figure_1_source_data.csv")
  saved <- paste0(out,c(".pdf",".png",".tiff"))
  if (any(!file.exists(saved)) || any(file.info(saved)$size <= 0)) stop("Figure export failed.",call. = FALSE)
  caption <- paste0(
    "Figure 1. Group-level amyloid-PET discrimination and individual-level classification agreement ",
    "across four plasma p-tau217 assays. A: frozen out-of-fold AUCs and 95% confidence intervals ",
    "using the same rank-AUC estimator and bootstrap as Table 2 in PET-matched anchor participants",
    if (length(pet_n) == 1L) paste0(" (n = ",pet_n,")") else "",
    ". B: each point represents one participant's smallest (horizontal) and largest (vertical) ",
    "Youden decision margins across four assays at the same blood draw. The pale red region ",
    "marks ranges crossing zero; red points indicate discordant negative and positive calls ",
    "(",k,"/",n,"; ",sprintf("%.1f",100*k/n),"%). Zero counts as positive. ",
    "The dashed diagonal denotes equal smallest and largest margins. C: four-assay binary ",
    "disagreement under triage, confirmatory and Youden rules; dots and horizontal bars show ",
    "estimates and full-pipeline participant-bootstrap 95% CIs, including calibration and threshold refitting. ",
    "D: mutually exclusive three-level classification patterns in complete-panel anchors; bars show ",
    "participant-bootstrap 95% CIs conditional on frozen classifications. Adjacent-only disagreement ",
    "excludes participants with any negative-positive pair. Negative-positive conflict occurred in ",
    conflict_n,"/",n," participants (",sprintf("%.1f",conflict$percent),"%). ",
    "Binary disagreement estimates and full-pipeline 95% CIs: ",
    binary_caption,". AUC describes discrimination against PET, not interchangeability ",
    "between the assays. The variance decomposition is reported in Figure S5.")
  writeLines(caption,file.path(results,"Figure_1_caption.txt"))
  if (interactive()) print(fig)
  message("Figure 1 saved: ",out,".{pdf,png,tiff}; all analysis inputs were computed in this run.")
})
