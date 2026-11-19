STUDY_B_ROOT_MANUAL <- Sys.getenv("STUDY_B_ROOT", unset = "")
if (nzchar(STUDY_B_ROOT_MANUAL) &&
    !file.exists(file.path(STUDY_B_ROOT_MANUAL, "CODE", "01_config.R")))
  stop("STUDY_B_ROOT points to an invalid project: ", STUDY_B_ROOT_MANUAL,
       call. = FALSE)
sb_require_manual <- function(objects, prior_scripts) {
  missing <- objects[!vapply(objects, exists, logical(1), inherits = TRUE)]
  if (length(missing)) {
    stop("Missing object(s): ", paste(missing, collapse = ", "),
         "\nRun ", prior_scripts,
         " manually first, in the same R session (Run sections or source the single script).",
         call. = FALSE)
  }
  invisible(TRUE)
}
sb_locate_root <- function() {
  frames <- sys.frames()
  ofiles <- vapply(frames, function(x) {
    if (is.null(x$ofile)) "" else as.character(x$ofile)[1]
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles)) {
    config_path <- normalizePath(tail(ofiles, 1), winslash = "/", mustWork = FALSE)
    code_dir <- dirname(config_path)
    if (basename(config_path) == "01_config.R" &&
        file.exists(file.path(code_dir, "01_config.R"))) return(dirname(code_dir))
  }
  env <- Sys.getenv("STUDY_B_ROOT", unset = "")
  if (nzchar(env) && file.exists(file.path(env, "CODE", "01_config.R"))) return(env)
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in seq_len(6)) {
    if (file.exists(file.path(d, "CODE", "01_config.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate the Study B package root.\n",
       "Source CODE/01_config.R or set STUDY_B_ROOT explicitly.",
       call. = FALSE)
}
sb_locate_data_dir <- function(root) {
  probe <- "UCBERKELEY_AMY_6MM.rda"
  env <- Sys.getenv("STUDY_B_DATA_DIR", unset = "")
  if (nzchar(env) && file.exists(file.path(env, probe))) return(env)
  candidates <- c(file.path(root, "01_DATA", "raw_rda"),
                  file.path(root, "DATA"),
                  file.path(root, "01_DATA"))
  for (cand in candidates) {
    if (file.exists(file.path(cand, probe))) return(cand)
  }
  candidates[1]
}
study_b_package_root <- sb_locate_root()
Sys.setenv(STUDY_B_ROOT = study_b_package_root)
if (!nzchar(Sys.getenv("STUDY_B_FIG3_REUSE_OBJECTS", unset = "")))
  Sys.setenv(STUDY_B_FIG3_REUSE_OBJECTS = "1")
study_b_code_dir <- file.path(study_b_package_root, "CODE")
study_b_data_dir <- sb_locate_data_dir(study_b_package_root)
study_b_dir      <- file.path(study_b_package_root, "RESULTS")
required_files <- c(
  file.path(study_b_data_dir, "FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.rda"),
  file.path(study_b_data_dir, "UCBERKELEY_AMY_6MM.rda")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing Study B input file(s):\n", paste(missing_files, collapse = "\n"),
       "\nSet STUDY_B_DATA_DIR to the folder containing the required .rda files, ",
       "or STUDY_B_ROOT/open STUDY_B_STANDALONE.Rproj first.", call. = FALSE)
}
dir.create(study_b_dir, recursive = TRUE, showWarnings = FALSE)
required_packages <- c(dplyr = "1.1.1", tidyr = "1.2.0", purrr = "1.0.0",
                       lme4 = "1.1-30", ggplot2 = "4.0.0",
                       lubridate = "1.9.0", scales = "1.2.0",
                       patchwork = "1.1.0")
sb_check_packages <- function(reqs) {
  bad <- character(0)
  for (p in names(reqs)) {
    if (!requireNamespace(p, quietly = TRUE)) {
      bad <- c(bad, sprintf("%s (not installed)", p))
    } else if (utils::packageVersion(p) < reqs[[p]]) {
      bad <- c(bad, sprintf("%s (have %s, need >= %s)", p,
                            as.character(utils::packageVersion(p)), reqs[[p]]))
    }
  }
  bad
}
bad_packages <- sb_check_packages(required_packages)
if (length(bad_packages)) {
  stop("Package requirements not met:\n  ", paste(bad_packages, collapse = "\n  "),
       "\nInstall/update the listed packages with install.packages() before running 01.", call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(lme4);  library(ggplot2); library(lubridate)
})
study_b_assay_map <- c(
  AlzPath   = "AlzPath_plasma_ptau217",
  C2N_ratio = "C2N_plasma_ptau217_ratio",
  Fujirebio = "Fuji_plasma_ptau217",
  Janssen   = "Janssen_plasma_ptau217"
)
study_b_assays <- names(study_b_assay_map)
study_b_seed              <- 20260903L
study_b_folds             <- 10L
study_b_centiloid_cutoff  <- 20
study_b_amyloid_window    <- 180
study_b_pet_qc_allowed    <- 2L
sb_seed <- function(label) {
  stopifnot(is.character(label), length(label) == 1L)
  b <- charToRaw(label)
  h <- sum(as.integer(b) * seq_along(b) * 7919)
  as.integer(study_b_seed + (h %% 100000L))
}
study_b_seed_folds  <- study_b_seed
study_b_seed_aarr   <- sb_seed("aarr")
study_b_seed_pair   <- sb_seed("pairwise")
study_b_seed_gray   <- sb_seed("gray_zone")
study_b_seed_spread <- sb_seed("assay_spread")
study_b_seed_dev    <- sb_seed("assay_deviation")
study_b_seed_vc     <- sb_seed("variance_components")
study_b_named_seeds <- c(
  folds = study_b_seed_folds, aarr = study_b_seed_aarr,
  pairwise = study_b_seed_pair, gray_zone = study_b_seed_gray,
  assay_spread = study_b_seed_spread, assay_deviation = study_b_seed_dev,
  variance_components = study_b_seed_vc)
if (anyDuplicated(study_b_named_seeds)) {
  stop("sb_seed() collision among named seeds: ",
       paste(names(study_b_named_seeds)[duplicated(study_b_named_seeds) |
                                          duplicated(study_b_named_seeds, fromLast = TRUE)],
             collapse = ", "),
       ". Rename one of the labels.", call. = FALSE)
}
study_b_boot_n       <- 2000L
study_b_boot_n_model <- 500L
study_b_boot_n_vc    <- 500L
study_b_boot_n_pipeline <- 500L
study_b_anchor_requires_pet <- FALSE
sb_make_core_contract <- function() {
  source_paths <- c(
    file.path(study_b_code_dir, c("01_config.R", "02_preprocess.R",
                                  "03_core_computation.R")),
    file.path(study_b_data_dir, c("FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.rda",
                                  "UCBERKELEY_AMY_6MM.rda"))
  )
  missing <- source_paths[!file.exists(source_paths)]
  if (length(missing)) {
    stop("Cannot construct the core-data contract; missing file(s):\n",
         paste(missing, collapse = "\n"), call. = FALSE)
  }
  list(
    contract_version = "study_b_core_v1",
    centiloid_cutoff = as.numeric(study_b_centiloid_cutoff),
    amyloid_window_days = as.integer(study_b_amyloid_window),
    pet_qc_allowed = paste(sort(unique(as.integer(study_b_pet_qc_allowed))),
                           collapse = ","),
    anchor_requires_pet = isTRUE(study_b_anchor_requires_pet),
    folds = as.integer(study_b_folds),
    seed = as.integer(study_b_seed),
    assay_map = paste(names(study_b_assay_map), unname(study_b_assay_map),
                      sep = "=", collapse = ";"),
    source_md5 = paste(basename(source_paths),
                       unname(tools::md5sum(source_paths)),
                       sep = "=", collapse = ";")
  )
}
study_b_sens_windows  <- c(180, 365)
study_b_sens_cutoffs  <- c(15, 20, 25, 30)
stopifnot(study_b_folds >= 2, study_b_boot_n >= 1,
          study_b_boot_n_pipeline >= 0, study_b_boot_n_vc >= 0,
          is.finite(study_b_centiloid_cutoff), study_b_amyloid_window > 0,
          is.logical(study_b_anchor_requires_pet),
          length(study_b_pet_qc_allowed) >= 1,
          all(is.finite(study_b_pet_qc_allowed)))
sb_load_rda <- function(path, expected_name = NULL) {
  e <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = e)
  if (!length(loaded_names)) stop("No object found in: ", path, call. = FALSE)
  if (!is.null(expected_name) && expected_name %in% loaded_names) {
    return(e[[expected_name]])
  }
  if (length(loaded_names) != 1L) {
    stop("Expected one object in ", path, "; found: ",
         paste(loaded_names, collapse = ", "),
         ". Supply expected_name explicitly.", call. = FALSE)
  }
  e[[loaded_names]]
}
sb_safe_num <- function(x, sentinel = TRUE) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x <- trimws(x)
    x[x %in% c("", "NA", "N/A", ".")] <- NA
  }
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  if (sentinel) x[x %in% c(-1, -4)] <- NA_real_
  x
}
sb_safe_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  d <- suppressWarnings(lubridate::ymd(x))
  bad <- is.na(d) & !is.na(x)
  if (any(bad)) d[bad] <- suppressWarnings(lubridate::mdy(x[bad]))
  d
}
sb_viskey <- function(v1, v2 = NULL) {
  a <- tolower(trimws(as.character(v1))); a[a == ""] <- NA_character_
  if (is.null(v2)) return(a)
  b <- tolower(trimws(as.character(v2))); b[b == ""] <- NA_character_
  dplyr::coalesce(b, a)
}
sb_pet_qc_pass <- function(qc_flag, allowed = study_b_pet_qc_allowed) {
  v <- suppressWarnings(as.numeric(as.character(qc_flag)))
  !is.na(v) & v %in% allowed
}
sb_qc_meaning <- function(qc_flag_raw, allowed = study_b_pet_qc_allowed) {
  keep <- ifelse(qc_flag_raw %in% allowed, "retained in primary",
                 "excluded from primary")
  base <- dplyr::case_when(
    qc_flag_raw ==  0 ~ "visual QC fail",
    qc_flag_raw ==  1 ~ "partial pass",
    qc_flag_raw ==  2 ~ "pass",
    qc_flag_raw == -1 ~ "not visually assessed",
    TRUE              ~ "other / undocumented")
  paste0(base, " (", keep, ")")
}
sb_first_nonmissing <- function(x) {
  z <- x[!is.na(x)]
  if (length(z)) z[1] else x[NA_integer_][1]
}
sb_dedup_visit <- function(dat, value_cols) {
  dat %>%
    mutate(.nonmiss = rowSums(!is.na(across(all_of(value_cols))))) %>%
    arrange(RID, date, desc(.nonmiss)) %>%
    distinct(RID, date, .keep_all = TRUE) %>%
    select(-.nonmiss)
}
sb_write <- function(x, file) {
  write.csv(x, file.path(study_b_dir, file), row.names = FALSE, na = "",
            fileEncoding = "UTF-8")
  invisible(x)
}
sb_have <- function(...) {
  nm <- c(...)
  all(vapply(nm, function(n) exists(n, envir = globalenv(), inherits = TRUE),
             logical(1)))
}
sb_group_folds <- function(ids, strata = NULL, k = study_b_folds,
                           seed = study_b_seed_folds) {
  if (is.null(strata)) strata <- rep("all", length(ids))
  stopifnot(length(ids) == length(strata), length(ids) > 0L,
            !anyNA(ids), length(k) == 1L, is.finite(k), k >= 2L, k == as.integer(k))
  if (any(tapply(as.character(strata), ids, function(z) length(unique(z))) > 1L))
    stop("Conflicting fold strata for the same participant.", call. = FALSE)
  strata[is.na(strata) | !nzchar(as.character(strata))] <- "unmatched"
  map <- tibble(RID = ids, stratum = as.character(strata)) %>%
    distinct(RID, .keep_all = TRUE)
  set.seed(seed)
  map %>%
    group_by(stratum) %>%
    slice_sample(prop = 1) %>%
    mutate(fold = rep(seq_len(k), length.out = n())) %>%
    ungroup() %>%
    select(RID, fold)
}
sb_boot_index <- function(data, id_col = "RID") {
  if (!id_col %in% names(data) || !nrow(data) || anyNA(data[[id_col]]))
    stop("Bootstrap requires nonempty data and nonmissing cluster IDs.", call. = FALSE)
  ids <- unique(data[[id_col]])
  idx <- split(seq_len(nrow(data)), factor(data[[id_col]], levels = ids))
  list(ids = ids, idx = idx, sizes = lengths(idx))
}
sb_boot_rows <- function(data, index, draw) {
  rows <- unlist(index$idx[draw], use.names = FALSE)
  d <- data[rows, , drop = FALSE]
  d$.boot_RID <- rep(seq_along(draw), index$sizes[draw])
  d
}
sb_cluster_boot <- function(data, statistic, B = study_b_boot_n,
                            seed = study_b_seed, id_col = "RID") {
  index <- sb_boot_index(data, id_col)
  n_id <- length(index$ids)
  set.seed(seed)
  out <- rep(NA_real_, B)
  if (B < 1) return(out[0])
  for (b in seq_len(B)) {
    draw <- sample.int(n_id, n_id, replace = TRUE)
    d <- sb_boot_rows(data, index, draw)
    out[b] <- tryCatch(as.numeric(statistic(d))[1], error = function(e) NA_real_)
  }
  out[is.finite(out)]
}
sb_cluster_boot_multi <- function(data, statistics, B = study_b_boot_n,
                                  seed = study_b_seed, id_col = "RID") {
  stopifnot(is.function(statistics))
  index <- sb_boot_index(data, id_col)
  n_id <- length(index$ids)
  probe <- statistics(data)
  k <- length(probe)
  out <- matrix(NA_real_, nrow = B, ncol = k,
                dimnames = list(NULL, names(probe)))
  if (B < 1) return(out[0, , drop = FALSE])
  set.seed(seed)
  for (b in seq_len(B)) {
    draw <- sample.int(n_id, n_id, replace = TRUE)
    d <- sb_boot_rows(data, index, draw)
    out[b, ] <- tryCatch(as.numeric(statistics(d)),
                         error = function(e) rep(NA_real_, k))
  }
  out
}
sb_ci <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(low = NA_real_, high = NA_real_))
  setNames(as.numeric(quantile(x, c(.025, .975), na.rm = TRUE)), c("low", "high"))
}
sb_boot_summary <- function(data, statistic, B, seed, id_col = "RID") {
  point <- tryCatch(as.numeric(statistic(data))[1], error = function(e) NA_real_)
  draws <- sb_cluster_boot(data, statistic, B = B, seed = seed, id_col = id_col)
  ci <- sb_ci(draws)
  c(estimate = point, low = unname(ci["low"]), high = unname(ci["high"]),
    bootstrap_successful = length(draws), bootstrap_attempted = B)
}
sb_agreement <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- as.integer(x[ok]); y <- as.integer(y[ok])
  if (!length(x)) return(c(agreement = NA_real_, kappa = NA_real_, gwet_ac1 = NA_real_))
  po <- mean(x == y); p1 <- mean(x == 1); p2 <- mean(y == 1)
  pe_k <- p1 * p2 + (1 - p1) * (1 - p2)
  pooled <- mean(c(x, y) == 1); pe_ac1 <- 2 * pooled * (1 - pooled)
  c(agreement = po,
    kappa    = ifelse(pe_k   == 1, NA_real_, (po - pe_k)   / (1 - pe_k)),
    gwet_ac1 = ifelse(pe_ac1 == 1, NA_real_, (po - pe_ac1) / (1 - pe_ac1)))
}
sb_pick_thresholds <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  if (any(!y[ok] %in% c(0, 1)))
    stop("Reference class must be coded 0/1.", call. = FALSE)
  y <- as.integer(y[ok]); score <- score[ok]
  if (length(unique(y)) < 2) stop("Training fold has one reference class", call. = FALSE)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) stop("Training fold has no positive or no negative", call. = FALSE)
  o <- order(score, decreasing = TRUE)
  s <- score[o]; yy <- y[o]
  tp <- cumsum(yy == 1); fp <- cumsum(yy == 0)
  last <- which(c(s[-1] != s[-length(s)], TRUE))
  cut         <- s[last]
  sensitivity <- tp[last] / n1
  specificity <- 1 - fp[last] / n0
  pick <- function(keep, ord) {
    if (!any(keep)) return(NA_real_)
    i <- which(keep)
    i[order(ord[[1]][i], ord[[2]][i], ord[[3]][i])][1]
  }
  zero <- rep(0, length(cut))
  i_tri <- pick(sensitivity >= .90, list(-specificity, -cut, zero))
  i_con <- pick(specificity >= .90, list(-sensitivity,  cut, zero))
  J <- sensitivity + specificity - 1
  i_you <- pick(rep(TRUE, length(cut)), list(-J, -sensitivity, -specificity))
  th <- c(triage = if (is.na(i_tri)) NA_real_ else cut[i_tri],
          confirmatory = if (is.na(i_con)) NA_real_ else cut[i_con],
          youden = if (is.na(i_you)) NA_real_ else cut[i_you])
  if (anyNA(th)) {
    stop("sb_pick_thresholds: no cut satisfied criterion(s): ",
         paste(names(th)[is.na(th)], collapse = ", "), call. = FALSE)
  }
  if (th[["triage"]] > th[["confirmatory"]]) {
    stop("sb_pick_thresholds: triage threshold (", signif(th[["triage"]], 4),
         ") exceeds confirmatory threshold (", signif(th[["confirmatory"]], 4),
         "); the three-class rule would be inverted. Independent 90%-sensitivity and ",
         "90%-specificity thresholds can overlap in a highly discriminating sample. ",
         "No threshold is reordered or changed automatically; review this training fold.",
         call. = FALSE)
  }
  th
}
sb_auc_rank <- function(y, score) {
  ok <- !is.na(y) & is.finite(score); y <- y[ok]; score <- score[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
sb_safe_ratio <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(NA_real_)
  num / den
}
sb_rules <- c(triage = "triage_positive",
              confirmatory = "confirmatory_positive",
              youden = "youden_positive")
sb_rule_matrix <- function(d, col, id_col = "RID") {
  keys <- d[, c(id_col, "assay")]
  if (anyDuplicated(keys)) {
    stop("sb_rule_matrix: ", sum(duplicated(keys)), " duplicated ", id_col,
         " x assay rows. Pass one observation per participant per assay ",
         "(e.g. sb_anchor), not repeated visits.", call. = FALSE)
  }
  w <- d %>%
    select(id = all_of(id_col), assay, value = all_of(col)) %>%
    pivot_wider(names_from = assay, values_from = value)
  missing_assay <- setdiff(study_b_assays, names(w))
  if (length(missing_assay)) {
    stop("sb_rule_matrix: missing assay column(s): ",
         paste(missing_assay, collapse = ", "), call. = FALSE)
  }
  as.matrix(w[, study_b_assays])
}
sb_complete_panel <- function(x) x[rowSums(is.na(x)) == 0, , drop = FALSE]
sb_discordant_rate <- function(x) {
  x <- sb_complete_panel(x)
  if (!nrow(x)) return(NA_real_)
  mean(apply(x, 1, function(z) length(unique(z)) > 1))
}
sb_aarr_stat <- function(col) function(d) {
  id_col <- if (".boot_RID" %in% names(d)) ".boot_RID" else "RID"
  sb_discordant_rate(sb_rule_matrix(d, col, id_col))
}
sb_pairs <- combn(study_b_assays, 2, simplify = FALSE)
sb_pair_stat <- function(col, p, which) function(d) {
  id_col <- if (".boot_RID" %in% names(d)) ".boot_RID" else "RID"
  dd <- filter(d, assay %in% p)
  keys <- dd[, c(id_col, "assay")]
  if (anyDuplicated(keys)) {
    stop("sb_pair_stat: duplicated ", id_col, " x assay rows.", call. = FALSE)
  }
  w <- dd %>%
    select(id = all_of(id_col), assay, value = all_of(col)) %>%
    pivot_wider(names_from = assay, values_from = value)
  unname(sb_agreement(w[[p[1]]], w[[p[2]]])[which])
}
cat("Study B standalone configuration loaded.\n",
    "Package: ", study_b_package_root, "\n",
    "Results: ", study_b_dir, "\n",
    "PET QC allowed (qc_flag): ", paste(study_b_pet_qc_allowed, collapse = ","), "\n",
    "Anchor requires PET: ", study_b_anchor_requires_pet, "\n", sep = "")
study_b_code_revision <- "2026-09-30-audit"
sb_code_label <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z[z %in% c("", "NA", "N/A", "UNKNOWN", "-1", "-4", "-5", "-9", "999")] <- NA_character_
  z
}
sb_map_diagnosis <- function(x) {
  z <- sb_code_label(x)
  dplyr::case_when(z %in% c("1", "CN", "CU", "NL", "NORMAL", "COGNITIVELY NORMAL") ~ "Cognitively unimpaired",
                   z %in% c("2", "MCI", "EMCI", "LMCI") ~ "MCI",
                   z %in% c("3", "AD", "DEMENTIA", "DEMENTIA DUE TO AD") ~ "Dementia",
                   TRUE ~ NA_character_)
}
sb_map_sex <- function(x) {
  z <- sb_code_label(x)
  dplyr::case_when(z %in% c("2", "F", "FEMALE") ~ TRUE,
                   z %in% c("1", "M", "MALE") ~ FALSE, TRUE ~ NA)
}
sb_map_race <- function(x) {
  z <- sb_code_label(x)
  dplyr::case_when(z %in% c("5", "WHITE") ~ "White",
                   z %in% c("4", "BLACK OR AFRICAN AMERICAN") ~ "Black or African American",
                   z %in% c("2", "ASIAN") ~ "Asian",
                   z %in% c("1", "3", "6", "AMERICAN INDIAN OR ALASKA NATIVE",
                            "AM INDIAN/ALASKAN", "AMERICAN INDIAN/ALASKA NATIVE",
                            "AMERICAN INDIAN OR ALASKAN NATIVE",
                            "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", "HAWAIIAN/OTHER PI",
                            "MORE THAN ONE", "MORE THAN ONE RACE", "MULTIPLE RACES",
                            "OTHER", "OTHER OR MULTIPLE RACE") ~ "Other or multiple race",
                   TRUE ~ NA_character_)
}
sb_map_ethnicity <- function(x) {
  z <- sb_code_label(x)
  dplyr::case_when(z %in% c("1", "HISPANIC OR LATINO", "HISP/LATINO",
                            "HISPANIC/LATINO") ~ TRUE,
                   z %in% c("2", "NOT HISPANIC OR LATINO", "NOT HISP/LATINO",
                            "NOT HISPANIC/LATINO") ~ FALSE, TRUE ~ NA)
}
sb_declared_unknown <- function(field, x) {
  z <- toupper(trimws(as.character(x)))
  common <- is.na(z) | z %in% c("", "NA", "N/A", "UNKNOWN", "UNK",
                                "NOT REPORTED", "REFUSED", "DECLINED",
                                "-1", "-4", "-5", "-9", "9", "999")
  common |
    (field == "PTRACCAT" & z %in% c("7")) |
    (field == "PTETHCAT" & z %in% c("3"))
}
sb_assert_fold_separation <- function(ev) {
  id <- if ("original_RID" %in% names(ev)) ev$original_RID else ev$RID
  if (anyNA(id) || anyNA(ev$fold)) stop("Missing participant identity or fold.")
  nfold <- tapply(ev$fold, id, function(z) length(unique(z)))
  if (any(nfold != 1L)) stop("Cross-fit leakage: copies of one ORIGINAL RID span folds.")
  invisible(TRUE)
}
sb_assign_boot_folds <- function(d, seed) {
  if (!"original_RID" %in% names(d)) stop("Bootstrap original_RID was not retained.")
  strata <- d %>% arrange(original_RID, plasma_date) %>% group_by(original_RID) %>%
    summarise(stratum = {z <- amyloid_positive[!is.na(amyloid_positive)];
    if (length(z)) as.character(first(z)) else "unmatched"}, .groups="drop")
  fm <- sb_group_folds(strata$original_RID, strata$stratum, seed=seed)
  names(fm)[names(fm)=="RID"] <- "original_RID"
  d$fold <- NULL
  d <- left_join(d, fm, by="original_RID", relationship="many-to-one")
  sb_assert_fold_separation(d)
  d
}
sb_lmer_converged <- function(m) {
  msg <- m@optinfo$conv$lme4$messages
  msg <- msg[!grepl("boundary.*singular", msg, ignore.case=TRUE)]
  code <- m@optinfo$conv$opt
  (is.null(code) || all(code == 0)) && length(msg)==0 &&
    is.finite(as.numeric(logLik(m))) && all(is.finite(lme4::getME(m,"theta")))
}
sb_model_audit <- new.env(parent=emptyenv())
sb_model_audit$rows <- list()
sb_checked_lmer <- function(formula, data, REML=TRUE, control=NULL, ...) {
  if (is.null(control)) control <- lme4::lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
  m <- lme4::lmer(formula,data=data,REML=REML,control=control,na.action=na.fail,...)
  retried <- FALSE
  if (!sb_lmer_converged(m)) {
    retried <- TRUE
    m <- lme4::lmer(formula,data=data,REML=REML,na.action=na.fail,
                    control=lme4::lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e6,rhoend=1e-8)),...)
  }
  ok <- sb_lmer_converged(m)
  dropped <- attr(lme4::getME(m, "X"), "col.dropped")
  if (length(dropped))
    stop("Mixed-model fixed effects are rank deficient: ",
         paste(names(dropped), collapse = ", "), ". Review the model subset.", call. = FALSE)
  sb_model_audit$rows[[length(sb_model_audit$rows)+1L]] <- data.frame(
    formula=paste(deparse(formula),collapse=" "), observations=stats::nobs(m),
    input_rows=nrow(data),
    converged=ok, singular=lme4::isSingular(m), retried=retried,
    message=paste(m@optinfo$conv$lme4$messages,collapse="; "))
  if (!ok) stop("Mixed model did not converge after retry; replicate/result rejected.")
  m
}
sb_model_sample_counts <- function(m) {
  mf <- stats::model.frame(m)
  required <- c("RID", "occ_id")
  if (!all(required %in% names(mf)) || anyNA(mf[required]))
    stop("Fitted model lacks complete RID/occ_id identity fields.", call. = FALSE)
  n_used <- stats::nobs(m)
  if (n_used != nrow(mf))
    stop("Model nobs differs from model.frame rows.", call. = FALSE)
  c(participants = dplyr::n_distinct(mf$RID),
    occasions = nrow(unique(mf[c("RID", "occ_id")])),
    n_observations = n_used)
}
sb_vc_fixed_terms <- function(data, assay_as_fixed=TRUE) {
  f <- if (dplyr::n_distinct(data$baseline_amyloid_group)>1)
    "years_from_baseline * baseline_amyloid_group" else "years_from_baseline"
  if (dplyr::n_distinct(data$phase)>1) f <- c(f,"phase")
  if (assay_as_fixed) f <- c(f,"assay")
  f
}
sb_deviation_fit <- function(d) {
  f <- assay_deviation ~ assay * years_from_baseline
  if (dplyr::n_distinct(d$phase)>1) f <- update(f,.~.+factor(phase))
  stats::lm(f,data=d,na.action=na.fail)
}
sb_clean_id <- function(x) {
  z <- trimws(as.character(x))
  z[is.na(z) | toupper(z) %in% c("", "NA", "N/A", "NULL", ".")] <- NA_character_
  z
}
sb_first_finite <- function(x) {
  z <- x[is.finite(x)]
  if (length(z)) z[1L] else NA_real_
}
sb_runtime_result_settings <- function() {
  list(
    bootstrap_counts = c(study_b_boot_n, study_b_boot_n_model,
                         study_b_boot_n_vc, study_b_boot_n_pipeline),
    sensitivity_windows = study_b_sens_windows,
    sensitivity_cutoffs = study_b_sens_cutoffs,
    extension_bootstrap = Sys.getenv(c("STUDY_B_UPGRADE_BOOT", "STUDY_B_B2_BOOT"),
                                     unset = "default"),
    laboratory_units = c(
      creatinine = tolower(Sys.getenv("STUDY_B_CREATININE_UNIT", unset = "mg/dL")),
      hemoglobin = tolower(Sys.getenv("STUDY_B_HEMOGLOBIN_UNIT", unset = "g/dL"))
    )
  )
}
.sb_expected_analysis_text_md5 <- "0dd9e69fc013dc02c3eb8a3ed7178a78"
.sb_expected_analysis_code_md5 <- c(
  "00_RUN_ALL.R" = "290dcb3fecb98131e764c72959691ee3",
  "02_preprocess.R" = "e2d69c81a1dc4d2f75ffe25e63f2471f",
  "03_core_computation.R" = "180a473d258941fbdad8c8c4c1f1775d",
  "04_table1_participant_characteristics.R" = "c4a2afd33e68212f7d589f9ac15befe9",
  "05_tableS1_sample_characteristics.R" = "fded3d0b1a8ccdbdf98429b141ee37ec",
  "06_figS1_flowchart.R" = "60ef9b9628d5aafdc1d3554a8d620c6c",
  "07_fig1_interchangeability.R" = "03b48294b4c72355c5145e78c1ceac49",
  "08_table2_clinical_performance.R" = "55249e4375e045fcf74470ba50a3e938",
  "09_figS5_variance_decomposition.R" = "dd7cb9b69b1c91e63f3e136da2324886",
  "10_tableS4_person_assay_audit.R" = "d1716978d6aa7d88ac49d2f2065fa1c9",
  "11_tableS9_longitudinal_summary.R" = "9a1a64bc4cdfb44a4e98acd1a8e4d465",
  "12_figS6_longitudinal_persistence.R" = "8a6ddf4c92f4a13ad0ca9edced8e46ba",
  "13_clinical_switching_and_transport.R" = "2079880d879c35b4eee2b45396722b29",
  "14_A2_pet_performance_exact_ci.R" = "4f91fc88714787f856b82878b58c559f",
  "15_A3_small_change_sensitivity.R" = "9bcf3a2a3cc1eee22c68f7b3755bcb9c",
  "17_figure3_longitudinal_consequences.R" = "f29061707e9315e16522161eec7cd211",
  "18_table3_translation_summary.R" = "48da80d2ef5d17b5f6c7651c90cd566c",
  "19_figS2_gray_zone_patterns.R" = "b62f84e88a4212088b1da17b9ecf9e72",
  "20_figS3_threshold_distributions.R" = "b58607ab2decb8af4e3cf3d986345267",
  "21_figS4_sensitivity_window_cutoff.R" = "d5929b4479d093b7bd88ad826e609b04",
  "22_tableS2_pairwise_agreement_full.R" = "84e65f8938c54026938b5faf4f1d023f",
  "23_tableS4_variance_components_full.R" = "424c2b4bd35ad5cf27bb03ddc885bae4",
  "24_tableS5_fixed_thresholds.R" = "7f603a500d555cbba1ea03643943c8c6",
  "25_tableS3_gray_zone_patterns_full.R" = "62c08572df7648407c00b77ebe8c4871",
  "26_tableS6_sensitivity_full.R" = "d48620d0397e2e678231cde5d4de3e02",
  "28_figS8_A3_small_change_sensitivity.R" = "79d0deb73cd825f2f4638fc022a96bf0",
  "29_B2_sample_key_audit.R" = "946b75e74cd699223a2ac3dce24e46b7",
  "30_B2_biological_modifiers.R" = "1586603bfe773fc9fce91ac1ef301d6b"
)
.sb_frozen_code_md5 <- c(
  "00_RUN_ALL.R" = "290dcb3fecb98131e764c72959691ee3",
  "01_config.R" = "ba45b7e79b8fa4af38c4c0ba3cf1c83c",
  "02_preprocess.R" = "e2d69c81a1dc4d2f75ffe25e63f2471f",
  "03_core_computation.R" = "180a473d258941fbdad8c8c4c1f1775d",
  "04_table1_participant_characteristics.R" = "c4a2afd33e68212f7d589f9ac15befe9",
  "05_tableS1_sample_characteristics.R" = "fded3d0b1a8ccdbdf98429b141ee37ec",
  "06_figS1_flowchart.R" = "60ef9b9628d5aafdc1d3554a8d620c6c",
  "07_fig1_interchangeability.R" = "03b48294b4c72355c5145e78c1ceac49",
  "08_table2_clinical_performance.R" = "55249e4375e045fcf74470ba50a3e938",
  "09_figS5_variance_decomposition.R" = "db2ab01da4f2962e1660230a7ed57f11",
  "10_tableS4_person_assay_audit.R" = "6815ddfaad5e4ede1bf590a114cdccae",
  "11_tableS9_longitudinal_summary.R" = "9a1a64bc4cdfb44a4e98acd1a8e4d465",
  "12_figS6_longitudinal_persistence.R" = "8a6ddf4c92f4a13ad0ca9edced8e46ba",
  "13_clinical_switching_and_transport.R" = "2079880d879c35b4eee2b45396722b29",
  "14_A2_pet_performance_exact_ci.R" = "4f91fc88714787f856b82878b58c559f",
  "15_A3_small_change_sensitivity.R" = "9bcf3a2a3cc1eee22c68f7b3755bcb9c",
  "16_figure2_switching_consequences.R" = "c18ce60e7130fe1c3a65773a79b77378",
  "17_figure3_longitudinal_consequences.R" = "f29061707e9315e16522161eec7cd211",
  "18_table3_translation_summary.R" = "48da80d2ef5d17b5f6c7651c90cd566c",
  "19_figS2_gray_zone_patterns.R" = "b62f84e88a4212088b1da17b9ecf9e72",
  "20_figS3_threshold_distributions.R" = "3433c804cb05c1b7440929eed99735f1",
  "21_figS4_sensitivity_window_cutoff.R" = "d5929b4479d093b7bd88ad826e609b04",
  "22_tableS2_pairwise_agreement_full.R" = "84e65f8938c54026938b5faf4f1d023f",
  "23_tableS4_variance_components_full.R" = "424c2b4bd35ad5cf27bb03ddc885bae4",
  "24_tableS5_fixed_thresholds.R" = "947f30231cda04645717c4e741f78929",
  "25_tableS3_gray_zone_patterns_full.R" = "62c08572df7648407c00b77ebe8c4871",
  "26_tableS6_sensitivity_full.R" = "d48620d0397e2e678231cde5d4de3e02",
  "27_headline_synthesis.R" = "74f29ea3e9f6cd9a5aa606bbbf2b81ac",
  "28_figS8_A3_small_change_sensitivity.R" = "79d0deb73cd825f2f4638fc022a96bf0",
  "29_B2_sample_key_audit.R" = "946b75e74cd699223a2ac3dce24e46b7",
  "30_B2_biological_modifiers.R" = "1586603bfe773fc9fce91ac1ef301d6b",
  "31_final_freeze_checks.R" = "a36244122f2fad0d003e576bacc40263",
  "study_b_display.R" = "5384b3c4aaafe27588662546477261ad"
)
.sb_frozen_file_md5 <- c(
  "01_amyloid_duplicate_audit.csv" = "b16ea523ec78eab01ce26696e6ce40c1",
  "01_amyloid_prefilter_audit.csv" = "5d46133126a884a655010042e7964103",
  "01_amyloid_qc_audit.csv" = "2a2af7474c4178e927219c26c9c285f9",
  "01_plasma_duplicate_audit.csv" = "089f3f7b13756f0e57cf63afc6901f34",
  "02_batch_field_audit.csv" = "6db4f43dc22b2bf7c190815b0360c45b",
  "02_cohort_manifest.csv" = "b6260cf87b56546768b836b68e8dbcc4",
  "02_crossfit_decision_scores.csv" = "7a04f112608d4504b853e999b82113bf",
  "02_training_fold_audit.csv" = "24a9067c1671765daec0bce0d328d4b2",
  "AUDIT_QC_occasion_counts.csv" = "15fcdc0e74ca925ff98fa3345b6c9285",
  "EXCLUDED_STALE_SUPPLEMENT_SOURCES.csv" = "5e1cb1581a7c3627d0302bd36f70c729",
  "FINAL_DISPLAY_EXPORT_CHECKS.csv" = "e65537c3db0946491645e09528b390a2",
  "FINAL_FREEZE_A3_STRICT_VS_LEGACY.csv" = "4819cc1fd896559d69bb0d3f64ffd4c6",
  "FINAL_FREEZE_CHECKS.csv" = "ef2fdeeba3886a7d3dedc46fc7acca19",
  "FINAL_MIXED_MODEL_DIAGNOSTICS.csv" = "c74ce44b67437b59ced1b20d621884d6",
  "Figure_1A_AARR.csv" = "3e6e8fc0032de16b62e26206af198b31",
  "Figure_1A_AARR_pipeline_bootstrap_draws.csv" = "795e67c57267f646a89f98462388207b",
  "Figure_1A_AARR_pipeline_bootstrap_failures.csv" = "74d6386ac534ca3faf0de052f4b8bfc0",
  "Figure_1B_pairwise_agreement_all_rules.csv" = "cea92b088216ee6e8ac0e3d8e55ee892",
  "Figure_1C_three_class_composition.csv" = "7f1ac244e5d74fadae8cbce1a73864c8",
  "Figure_1_AUC.csv" = "ad741beb418171412913e3f4166ea49c",
  "Figure_1_interchangeability.pdf" = "e2c8b115dbc5a39583a2fe672d9a97f9",
  "Figure_1_interchangeability.png" = "5336e07d19ca6b40cce6ca9e3adcc102",
  "Figure_1_interchangeability.tiff" = "00506fa9d0af1aaf7c55eec48ce12a3f",
  "Figure_1_participant_panel.csv" = "ca0e26edfb130ab767fb252a3d04602a",
  "Figure_1_source_data.csv" = "49c5a7b926e108b5d8904b28e59a983f",
  "Figure_2_source_data.csv" = "1fa117c08951e1592c982f32ff25ae91",
  "Figure_2_switching_consequences.pdf" = "8439733e6238df87e484b67c2a1f1e73",
  "Figure_2_switching_consequences.png" = "72cc5edcdf07162390ac479774451733",
  "Figure_2_switching_consequences.tiff" = "857cf47e68835e98ecce8336fb99d632",
  "Figure_3_distribution_source_data.csv" = "fe6491f17450abf8724131831ac20955",
  "Figure_3_longitudinal_switch_consequences.pdf" = "4dd9e2205334f0efc867d759b68cf6ae",
  "Figure_3_longitudinal_switch_consequences.png" = "c355251435e302d16e77e571c066a47a",
  "Figure_3_longitudinal_switch_consequences.tiff" = "31a849172a863980711ed74055274b0f",
  "Figure_3_source_data.csv" = "f37357a59bdb5ce533e95a147038c320",
  "Figure_S1_PET_selection_audit.csv" = "12b72953520151aee1646c7f6bad060a",
  "Figure_S1_collection_disposition.csv" = "ce562dd0948175827f16ee38c515dfa0",
  "Figure_S1_flowchart.pdf" = "91c64a1434f8306e5c6bf93dbfc1eb2a",
  "Figure_S1_flowchart.png" = "cec66b3f36e88eb557278ae800d31f41",
  "Figure_S1_flowchart.tiff" = "7c696f06c426f82c126482f7965f63b6",
  "Figure_S1_flowchart_numbers.csv" = "18a9e27a84826d5fefa67e9743683a1a",
  "Figure_S1_plasma_selection_audit.csv" = "6f686a7cf22c1d27f6b1c2565dca3ec0",
  "Figure_S2_gray_zone_patterns_top15.pdf" = "9ba280571f8dd88b86a3d9d4251acf6f",
  "Figure_S2_gray_zone_patterns_top15.png" = "cb0025ae68562cf3f449656398b52686",
  "Figure_S2_gray_zone_patterns_top15.tiff" = "a589509a79ae5fdb8963bd1e4c1a2588",
  "Figure_S3_full_sample_fixed_thresholds.csv" = "92e537f3e2ebd9e888fb52e8f0eb9bc8",
  "Figure_S3_threshold_distributions.pdf" = "fc4bb960b59cd9907762a4aec05386a7",
  "Figure_S3_threshold_distributions.png" = "6335804b5460097b5e25da495394ca58",
  "Figure_S3_threshold_distributions.tiff" = "f843d2066332d22c3977fb88c67bc73b",
  "Figure_S4_sensitivity_ratio.pdf" = "2878c94082144373d567335287fe2f7d",
  "Figure_S4_sensitivity_ratio.png" = "91d3af43fa7c4e881b4141e7a1f7457a",
  "Figure_S4_sensitivity_ratio.tiff" = "a2f03fcdaa0d1fcc6ebded5ee7154649",
  "Figure_S4_sensitivity_window_cutoff.csv" = "94c497d8e1ef27451dd3c85e1b7a1e1d",
  "Figure_S5_assay_switch_flip_by_margin.csv" = "47db89623a3671fb5742ed7db187e5aa",
  "Figure_S5_model_fit_audit.csv" = "6048146569c4dac43e9ee7ac6c05521e",
  "Figure_S5_person_by_assay_conditional_modes.csv" = "fba9889b6ccb8321da8602965d036de6",
  "Figure_S5_variance_components_main.csv" = "b4d2e07a8745d5c174f4ff6590d7b7e0",
  "Figure_S5_variance_decomposition.pdf" = "5b2d2bd5ce83fdde199add5e962f5dc1",
  "Figure_S5_variance_decomposition.png" = "2c58d0255265920442312aab72206227",
  "Figure_S5_variance_decomposition.tiff" = "d87d9f490d07fe2741cd13f18aba23b2",
  "Figure_S5_variance_ratio_bootstrap_audit.csv" = "f0acc8380ad68f8d1375cb1ca3ed6814",
  "Figure_S5_variance_ratio_bootstrap_draws.csv" = "014c0caa8bb8be89c27cd8b5cedec8ed",
  "Figure_S5_variance_ratio_bootstrap_failures.csv" = "a2f664bcc450b498f1eecb7f261c6665",
  "Figure_S5_variance_ratio_ci.csv" = "0b1ee9d4e0f4ddb198c86243f2dcbe36",
  "Figure_S6_deviation_first_vs_last.csv" = "ce21459f593845c4f021de181d37cca2",
  "Figure_S6_deviation_first_vs_last.png" = "7cf129394630cb2a0638862dd0433709",
  "Figure_S6_deviation_persistence_by_assay.csv" = "731dd8e389caa950a637478c5ff98fd0",
  "Figure_S6_deviation_sign_transition.csv" = "5e5fdd357b7045a638f90bd7d5f1155c",
  "Figure_S6_longitudinal_persistence.pdf" = "935c2184714a492f64d8142946105aac",
  "Figure_S6_longitudinal_persistence.png" = "3caa8e868227b7c8b2603b34ecf723f1",
  "Figure_S6_longitudinal_persistence.tiff" = "76458cca5add483e39a3191df525ad03",
  "Figure_S7_remeasurement_vs_switching.csv" = "da461e41517726cbe17943d881c9a3e7",
  "Figure_S7_remeasurement_vs_switching.pdf" = "29b2879d66f3264009b1d3a7dab03b71",
  "Figure_S7_remeasurement_vs_switching.png" = "8e3f9069311b949416f8f19a5e188225",
  "Figure_S7_remeasurement_vs_switching.tiff" = "fd5d8174419f686409b590eecf32f7ed",
  "Figure_S8_A3_small_change_sensitivity.pdf" = "7dfd73c408b312f141179630d3e270ae",
  "Figure_S8_A3_small_change_sensitivity.png" = "1fa082ad9095e4aaaff874625127cedb",
  "Figure_S8_A3_small_change_sensitivity.tiff" = "f109cea9f57e7f8d3e9928a76bb5a11a",
  "Figure_S8_source_data.csv" = "33742681d501fa338b5e553cc72cc250",
  "HEADLINE_analysis_upgrades.csv" = "94683c465b0b72803ff3d61b6651c60e",
  "HEADLINE_results.csv" = "fc3e116364a9d7f8ed37f72fe6bd085c",
  "MANUSCRIPT_DISPLAY_MAP.csv" = "f921d5c860c4e3ebb983e761f09edadf",
  "SUPPLEMENT_TABLE_INDEX.csv" = "869186e29847450bcd70344e3d8bbe40",
  "Table_1_denominators.csv" = "1f763c3898fabb4d878861943843ec32",
  "Table_1_encoding_audit.csv" = "31b17179878aa27f5c5725bb59315dc9",
  "Table_1_participant_audit.csv" = "0a18af41f380748e967c189c88fd90aa",
  "Table_1_participant_characteristics.csv" = "df80ae71a6708daf7d7f069b185530da",
  "Table_1_source_audit.csv" = "6a3968d71d375eb83d2a64ee4b92e8b7",
  "Table_2_crossfit_clinical_performance.csv" = "523da7201f029c1e13df526bc095d60c",
  "Table_3_translation_consequences.csv" = "b8ee18048380639d40c73ac449ef9d45",
  "Table_S10_B1_PET_direction_binary_switch.csv" = "3ccfebb0399f519381d39267dafbede7",
  "Table_S10_B1_PET_direction_three_class_switch.csv" = "a03edd6e7393cc183486284acffae6ed",
  "Table_S10_B3_hypothetical_prevalence_by_assay.csv" = "31e43d7377345085d8fab95d8cb7db9a",
  "Table_S10_B3_hypothetical_prevalence_switching.csv" = "9c5dec625ad2bb55c99951c847cc9aed",
  "Table_S10_final.csv" = "2099ec04b20df960007dadd599b1ddbe",
  "Table_S11_B2_bootstrap_diagnostics.csv" = "0215567c3794fdf70352a36569cb45af",
  "Table_S11_B2_bootstrap_replicates.csv" = "986ba3c54e606a30f8da1d0e85786287",
  "Table_S11_B2_matched_measurement_audit.csv" = "94761ab9489f147a26b0c9b770a43689",
  "Table_S11_B2_modifier_coverage.csv" = "f2cb5230bda48f44a34d732c2b45ce85",
  "Table_S11_B2_modifier_value_audit.csv" = "caad7acae853ba527d6229a6503699b8",
  "Table_S11_B2_sample_key_audit.csv" = "43986b6d4b9f09a15e8f22a58ef2d282",
  "Table_S11_B2_unit_contract.csv" = "ad2864928a9f0ef6e7d680799d41a023",
  "Table_S11_B2_variance_explanation.csv" = "bb6a5158b2eb173a6a0ce1a4b88c68f3",
  "Table_S11_B2_vitals_unit_audit.csv" = "2a17fa381cfb89f623d167adcf2fbbc9",
  "Table_S11_final.csv" = "07aa80d602df85fd5bdc3ae1bf02a71b",
  "Table_S1_assay_decision_prevalence.csv" = "5fd8adafed3ef8edd805139a80b5a2f4",
  "Table_S1_final.csv" = "969a73c5fe65df12abffe1d378b478eb",
  "Table_S1_three_class_prevalence.csv" = "72c7919cd1c178964c089a63da8ba20e",
  "Table_S2_final.csv" = "9cfeaba137ae10ffe3f092dace0b9b0e",
  "Table_S2_pairwise_agreement_full.csv" = "cea92b088216ee6e8ac0e3d8e55ee892",
  "Table_S3_final.csv" = "ca5bede8d558798bba169ffac9ac753f",
  "Table_S3_gray_zone_pairwise.csv" = "58acf2a5527d31520fd0287a98e08ff7",
  "Table_S3_gray_zone_patterns_full.csv" = "ff3cdddafdd7953ca28119fabfe275ae",
  "Table_S3_pattern_coverage.csv" = "567b4c235deb8b6678746eb586b4498b",
  "Table_S4_components_full.csv" = "3c9e40d59f55fcb884f792e938b7f051",
  "Table_S4_final.csv" = "d8132d3fb818185026d3338d5c18cae5",
  "Table_S4_fixed_threshold_audit.csv" = "92e537f3e2ebd9e888fb52e8f0eb9bc8",
  "Table_S4_person_assay_audit.csv" = "cbea0615fba28596b7b5c6295427df68",
  "Table_S4_variance_components_full.csv" = "61c9054b59708460f1857fd0c5e6ad2e",
  "Table_S5_final.csv" = "159322e31a3e961f1fd240903d67af14",
  "Table_S5_fixed_thresholds.csv" = "92e537f3e2ebd9e888fb52e8f0eb9bc8",
  "Table_S6_final.csv" = "defcb232fe1c9aca95a773e719b9c23b",
  "Table_S6_sensitivity_window_cutoff_full.csv" = "94c497d8e1ef27451dd3c85e1b7a1e1d",
  "Table_S7_A1_assay_switching_consequences.csv" = "740e882e01c3ee0dd7b1850482d8d8c3",
  "Table_S7_A1_switching_risk_by_source_margin.csv" = "a67fba45d11970b6e9d4398ec2bc87b6",
  "Table_S7_A1_switching_risk_by_source_margin_CI.csv" = "e2136c7e676330f4e12f14e9d989d086",
  "Table_S7_final.csv" = "fe9f1ea00f9b8a3f9ba8698487747aca",
  "Table_S8_A2_PET_performance_exact_CI.csv" = "a1c3b992c0baab55df2001d39464cd4e",
  "Table_S8_A2_symptomatic_PET_performance.csv" = "e89cb714c2dc35b0b6be97ce8b41fd45",
  "Table_S8_A2_symptomatic_disagreement.csv" = "5895586805d145fd857c4c54c2a193e8",
  "Table_S8_A2_symptomatic_switching_consequences.csv" = "484efe24f08159627fb46184f521618a",
  "Table_S8_A2_symptomatic_variance_components.csv" = "c23bd45dbc0d661a481005a41e4eecb0",
  "Table_S8_final.csv" = "9768e012bddd38e328721ed78ed280ed",
  "Table_S9_A3_date_alignment_audit.csv" = "f0acd35e8a13c56c9d938dbd118af988",
  "Table_S9_A3_direction_reversal_zero_change_audit.csv" = "7053f60b17816b37ea583d5f48bb65ad",
  "Table_S9_A3_longitudinal_pseudo_switch.csv" = "26f22341e4ce5e2a65ee53bcb666fe85",
  "Table_S9_A3_longitudinal_pseudo_switch_magnitudes.csv" = "e965da6eef878ca61086974557a02749",
  "Table_S9_A3_same_change_quartile_sensitivity.csv" = "fe49ace1abce3cd2a8921b196afa8444",
  "Table_S9_A3_small_change_threshold_sensitivity.csv" = "e8c928c00fa18c6e326baa70fa04efb8",
  "Table_S9_final.csv" = "e314488a178221211c945b659c41b129",
  "Table_S9_supporting_deviation_by_assay.csv" = "07205c17d452313b67483ea82a519ab9",
  "Table_S9_supporting_longitudinal_persistence.csv" = "718bc749aff0e8938c13ee4d4c2fbbbe",
  "Table_S9_supporting_persistence_calibration_audit.csv" = "cec5ad5d04804b0b65a78211473a1a7e"
)
.sb_frozen_preprocess_md5 <- "8d48d67dd8f92b768bd7e22adbc2647a"
.sb_frozen_core_md5 <- "8847220e70d910ca4311c3299b7ee68b"
.sb_reporting_scripts <- c("01_config.R", "16_figure2_switching_consequences.R",
                          "27_headline_synthesis.R", "31_final_freeze_checks.R",
                          "study_b_display.R", "32_reporting_input_audit.R")
.sb_analysis_text_md5 <- function(path) {
  z <- readLines(path, encoding = "UTF-8", warn = FALSE)
  part <- function(start, end) {
    a <- which(startsWith(z, start)); b <- which(startsWith(z, end))
    if (length(a) != 1L || !any(b > a))
      stop("Cannot identify audited analysis block in 01_config.R.", call. = FALSE)
    b <- b[b > a][1L]
    z[seq.int(a, b - 1L)]
  }
  blocks <- c(part("study_b_assay_map <-", "sb_make_core_contract <-"),
              part("## Sensitivity grids used", "sb_model_audit <-"),
              part("sb_checked_lmer <-", "  sb_model_audit$rows[["),
              part("sb_vc_fixed_terms <-", "sb_runtime_result_settings <-"))
  p <- tempfile(); on.exit(unlink(p), add = TRUE)
  con <- file(p, "wb"); writeLines(enc2utf8(blocks), con, useBytes = TRUE); close(con)
  unname(tools::md5sum(p))
}
.sb_source_hashes <- function(x) {
  p <- strsplit(x, ";", fixed = TRUE)[[1]]
  out <- sub("^[^=]+=", "", p)
  names(out) <- sub("=.*$", "", p)
  out
}
.sb_core_equivalent <- function(old, now) {
  if (!is.list(old) || !is.list(now)) return(FALSE)
  a <- old; b <- now; a$source_md5 <- b$source_md5 <- NULL
  if (!identical(a, b)) return(FALSE)
  h_old <- .sb_source_hashes(old$source_md5)
  h_now <- .sb_source_hashes(now$source_md5)
  keep <- setdiff(names(h_old), "01_config.R")
  identical(sort(names(h_old)), sort(names(h_now))) &&
    identical(unname(h_old[keep]), unname(h_now[keep])) &&
    identical(unname(h_old["01_config.R"]), .sb_frozen_code_md5[["01_config.R"]])
}
.sb_check_audited_analysis <- function() {
  if (!identical(.sb_analysis_text_md5(file.path(study_b_code_dir, "01_config.R")),
                 .sb_expected_analysis_text_md5))
    stop("Analysis helpers/settings in 01 differ from the audited version. Frozen-result compatibility rejected.", call. = FALSE)
  paths <- file.path(study_b_code_dir, names(.sb_expected_analysis_code_md5))
  if (!all(file.exists(paths)) ||
      !identical(unname(tools::md5sum(paths)), unname(.sb_expected_analysis_code_md5)))
    stop("An analysis script differs from the audited version. Regenerate its affected results; reporting-only reuse is rejected.", call. = FALSE)
  invisible(TRUE)
}
.sb_check_audited_analysis()
sb_make_core_contract_native <- sb_make_core_contract
.sb_active_frozen_core <- NULL
.anchor_meta <- file.path(study_b_dir, "Table_S7_A1_assay_switching_consequences.csv.provenance.rds")
.anchor_csv <- sub("[.]provenance[.]rds$", "", .anchor_meta)
if (file.exists(.anchor_meta) && file.exists(.anchor_csv) &&
    identical(unname(tools::md5sum(.anchor_csv)),
              .sb_frozen_file_md5[[basename(.anchor_csv)]])) {
  .origin <- readRDS(.anchor_meta)
  if (identical(.origin$contract$code_md5, .sb_frozen_code_md5) &&
      .sb_core_equivalent(.origin$contract$core, sb_make_core_contract_native()))
    .sb_active_frozen_core <- .origin$contract$core
}
sb_make_core_contract <- function() {
  now <- sb_make_core_contract_native()
  if (!is.null(.sb_active_frozen_core) && .sb_core_equivalent(.sb_active_frozen_core, now))
    return(.sb_active_frozen_core)
  now
}
.sb_current_code <- function() {
  h <- tools::md5sum(sort(list.files(study_b_code_dir, pattern = "[.]R$", full.names = TRUE)))
  names(h) <- basename(names(h)); h
}
sb_result_contract_used <- list(
  provenance_version = "study_b_reporting_v2",
  core = sb_make_core_contract(),
  code_md5 = .sb_current_code(),
  analysis_text_md5 = .sb_expected_analysis_text_md5,
  analysis_code_md5 = .sb_expected_analysis_code_md5,
  settings = sb_runtime_result_settings()
)
sb_assert_code_current <- function() {
  if (!identical(.sb_current_code(), sb_result_contract_used$code_md5))
    stop("CODE changed during this R session. Run the revised 01 again after all replacements; do not rerun statistical producers for a reporting-only change.", call. = FALSE)
  .sb_check_audited_analysis()
  if (!identical(sb_runtime_result_settings(), sb_result_contract_used$settings) ||
      !identical(sb_make_core_contract(), sb_result_contract_used$core))
    stop("Analysis settings or raw/core inputs changed. Reporting-only reuse rejected.", call. = FALSE)
  invisible(TRUE)
}
.sb_reporting_state <- new.env(parent = emptyenv())
.sb_reporting_state$producer <- NULL
.sb_reporting_state$inputs <- character()
sb_begin_reporting_step <- function(script) {
  sb_assert_code_current()
  if (!script %in% .sb_reporting_scripts || script == "01_config.R")
    stop("Unregistered reporting step: ", script, call. = FALSE)
  .sb_reporting_state$producer <- script
  .sb_reporting_state$inputs <- character()
  invisible(TRUE)
}
.sb_resolve_input <- function(path) {
  if (file.exists(path)) return(path)
  candidates <- c(file.path(study_b_data_dir, basename(path)),
                  file.path(study_b_dir, basename(path)))
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(hit[1])
  path
}
.sb_verify_hashes <- function(h) {
  if (!length(h)) return(TRUE)
  if (is.null(names(h))) return(FALSE)
  p <- vapply(names(h), .sb_resolve_input, character(1))
  all(file.exists(p)) && identical(unname(tools::md5sum(p)), unname(h))
}
.sb_note_input <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  .sb_reporting_state$inputs[path] <- unname(tools::md5sum(path))
  invisible(path)
}
.sb_settings_match <- function(old) {
  if (!is.list(old)) return(FALSE)
  if (is.null(old$laboratory_units))
    old$laboratory_units <- c(creatinine = "mg/dl", hemoglobin = "g/dl")
  identical(old, sb_runtime_result_settings())
}
sb_validate_result <- function(file, seen = character()) {
  sb_assert_code_current()
  if (identical(basename(file), file)) file <- file.path(study_b_dir, file)
  if (!file.exists(file) || !file.exists(paste0(file, ".provenance.rds")))
    stop("Result or provenance absent: ", file, call. = FALSE)
  path <- normalizePath(file, winslash = "/", mustWork = TRUE)
  if (path %in% seen) stop("Circular result provenance: ", file, call. = FALSE)
  meta <- readRDS(paste0(file, ".provenance.rds"))
  hash <- unname(tools::md5sum(file))
  if (!identical(hash, meta$result_md5) ||
      !identical(meta$contract$core, sb_result_contract_used$core) ||
      !.sb_settings_match(meta$contract$settings) ||
      !.sb_verify_hashes(meta$extra_input_md5))
    stop("Changed/stale result, source data or settings: ", basename(file), call. = FALSE)
  if (identical(meta$provenance_version, "study_b_reporting_v2")) {
    producer <- meta$producer_script
    if (!is.character(producer) || length(producer) != 1L ||
        !producer %in% .sb_reporting_scripts ||
        !identical(meta$producer_code_md5, unname(tools::md5sum(file.path(study_b_code_dir, producer)))) ||
        !identical(meta$contract$analysis_text_md5, .sb_expected_analysis_text_md5) ||
        !identical(meta$contract$analysis_code_md5, .sb_expected_analysis_code_md5) ||
        !.sb_verify_hashes(meta$input_md5))
      stop("Reporting output or its producer/input changed: ", basename(file), call. = FALSE)
    if (!is.null(meta$display_helper_md5) &&
        !identical(meta$display_helper_md5, unname(tools::md5sum(file.path(study_b_code_dir, "study_b_display.R")))))
      stop("Figure display helper changed; redraw: ", basename(file), call. = FALSE)
    for (p in names(meta$input_md5)) {
      p <- .sb_resolve_input(p)
      if (file.exists(paste0(p, ".provenance.rds")))
        sb_validate_result(p, seen = c(seen, path))
    }
  } else if (identical(meta$provenance_version, "study_b_computation_v2")) {
    if (!identical(meta$contract, sb_result_contract_used))
      stop("New computational result has a different execution contract: ", basename(file), call. = FALSE)
  } else {
    expected <- .sb_frozen_file_md5[basename(file)]
    if (length(expected) != 1L || is.na(expected) ||
        !identical(hash, unname(expected)) ||
        !identical(meta$contract$code_md5, .sb_frozen_code_md5))
      stop("This legacy result is not in the audited frozen release: ", basename(file), call. = FALSE)
  }
  invisible(meta)
}
sb_stamp_result <- function(path, extra_inputs = character()) {
  sb_assert_code_current()
  if (is.null(.sb_reporting_state$producer)) {
    stopifnot(file.exists(path), all(file.exists(extra_inputs)))
    meta <- list(provenance_version = "study_b_computation_v2",
                 contract = sb_result_contract_used,
                 result_md5 = unname(tools::md5sum(path)),
                 extra_input_md5 = tools::md5sum(extra_inputs))
    saveRDS(meta, paste0(path, ".provenance.rds"))
    return(invisible(path))
  }
  stopifnot(file.exists(path), all(file.exists(extra_inputs)))
  input_hashes <- .sb_reporting_state$inputs
  if (length(input_hashes) && normalizePath(path, winslash = "/", mustWork = TRUE) %in% names(input_hashes))
    stop("Refusing to overwrite a frozen input consumed by this step: ", path, call. = FALSE)
  meta <- list(
    provenance_version = "study_b_reporting_v2",
    contract = sb_result_contract_used,
    producer_script = .sb_reporting_state$producer,
    producer_code_md5 = unname(tools::md5sum(file.path(study_b_code_dir, .sb_reporting_state$producer))),
    result_md5 = unname(tools::md5sum(path)),
    input_md5 = input_hashes,
    extra_input_md5 = tools::md5sum(extra_inputs)
  )
  if (tolower(tools::file_ext(path)) %in% c("pdf", "png", "tiff", "tif"))
    meta$display_helper_md5 <- unname(tools::md5sum(file.path(study_b_code_dir, "study_b_display.R")))
  saveRDS(meta, paste0(path, ".provenance.rds")); invisible(path)
}
sb_save_csv <- function(x, file, ...) {
  sb_assert_code_current()
  path <- normalizePath(file, winslash = "/", mustWork = FALSE)
  if (path %in% names(.sb_reporting_state$inputs))
    stop("Refusing to overwrite a frozen input: ", file, call. = FALSE)
  utils::write.csv(x, file = file, ...); sb_stamp_result(file); invisible(x)
}
sb_write <- function(x, file) {
  sb_save_csv(x, file.path(study_b_dir, file), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
sb_read_result <- function(file, ...) {
  if (identical(basename(file), file)) file <- file.path(study_b_dir, file)
  sb_validate_result(file); .sb_note_input(file); utils::read.csv(file, ...)
}
sb_read_core <- function(path = file.path(study_b_dir, "02_core_results.rds")) {
  sb_assert_code_current()
  x <- readRDS(path)
  need <- c("sb_scored", "sb_anchor", "sb_rep_dat", "sb_events", "sb_decision_long", "sb_core_contract_used")
  if (!all(need %in% names(x)) || !identical(x$sb_core_contract_used, sb_make_core_contract()))
    stop("Stale/incomplete core snapshot.", call. = FALSE)
  if (!is.null(.sb_active_frozen_core) &&
      !identical(unname(tools::md5sum(path)), .sb_frozen_core_md5))
    stop("Frozen core snapshot differs from the audited snapshot.", call. = FALSE)
  .sb_note_input(path); x
}
sb_end_reporting_step <- function() {
  .sb_reporting_state$producer <- NULL
  .sb_reporting_state$inputs <- character()
  invisible(TRUE)
}
.display_helper <- file.path(study_b_code_dir,"study_b_display.R")
if (!file.exists(.display_helper)) stop("Missing CODE/study_b_display.R; replace the complete CODE folder.",call. = FALSE)
source(.display_helper,local = .GlobalEnv,encoding = "UTF-8",chdir = FALSE)
