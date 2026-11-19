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
.table1_config_ok <- function() {
  required <- c(
    "sb_have", "sb_load_rda", "sb_safe_date", "sb_safe_num",
    "sb_first_nonmissing", "sb_code_label", "sb_map_sex", "sb_map_race",
    "sb_map_ethnicity", "sb_map_diagnosis", "sb_declared_unknown", "sb_write",
    "sb_make_core_contract"
  )
  if (!all(vapply(required, exists, logical(1), mode = "function",
                  inherits = TRUE))) return(FALSE)
  identical(
    unname(sb_map_race("American Indian or Alaskan Native")),
    "Other or multiple race"
  )
}
if (!.table1_config_ok()) {
  stop("Run 01_config.R manually first, in this same R session.", call. = FALSE)
}
if (!.table1_config_ok()) {
  stop(
    "03_table1 requires a compatible 01_config.R. Replace both supplied files ",
    "together; the race crosswalk must recognize ",
    "'American Indian or Alaskan Native'.",
    call. = FALSE
  )
}
.table1_core_path <- file.path(study_b_dir, "02_core_results.rds")
if (!file.exists(.table1_core_path)) {
  stop("Table 1 requires RESULTS/02_core_results.rds. Run 03_core_computation.R ",
       "Run scripts 01-03 manually first.", call. = FALSE)
}
.table1_core <- readRDS(.table1_core_path)
.table1_core_required <- c("sb_anchor", "sb_events", "cohort",
                           "sb_core_contract_used")
.table1_core_missing <- setdiff(.table1_core_required, names(.table1_core))
if (length(.table1_core_missing)) {
  stop("The core snapshot is incomplete or predates the consistency fix; missing: ",
       paste(.table1_core_missing, collapse = ", "),
       ". Rerun scripts 01-04 manually from a clean session.", call. = FALSE)
}
if (!identical(.table1_core$sb_core_contract_used, sb_make_core_contract())) {
  stop("Table 1 refused a stale core snapshot: configuration, upstream code, ",
       "or input files have changed since 02_core_results.rds was created. ",
       "Rerun scripts 01-04 manually.", call. = FALSE)
}
list2env(.table1_core[.table1_core_required], envir = .GlobalEnv)
table1_files <- c(
  PTDEMOG = file.path(study_b_data_dir, "PTDEMOG.rda"),
  APOERES = file.path(study_b_data_dir, "APOERES.rda"),
  DXSUM   = file.path(study_b_data_dir, "DXSUM.rda")
)
missing_table1_files <- table1_files[!file.exists(table1_files)]
if (length(missing_table1_files)) {
  stop("Table 1 cannot be generated. Missing required source file(s):\n",
       paste(missing_table1_files, collapse = "\n"), call. = FALSE)
}
ptdemog <- sb_load_rda(table1_files[["PTDEMOG"]], "PTDEMOG")
apoeres <- sb_load_rda(table1_files[["APOERES"]], "APOERES")
dxsum   <- sb_load_rda(table1_files[["DXSUM"]],   "DXSUM")
required_columns <- list(
  PTDEMOG = c("RID", "VISDATE", "PTGENDER", "PTDOB", "PTEDUCAT",
              "PTRACCAT", "PTETHCAT"),
  APOERES = c("RID", "GENOTYPE"),
  DXSUM   = c("RID", "EXAMDATE", "DIAGNOSIS")
)
for (nm in names(required_columns)) {
  dat <- get(tolower(nm))
  miss <- setdiff(required_columns[[nm]], names(dat))
  if (length(miss)) {
    stop(nm, " is missing required field(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
}
.table1_find_data_file <- function(stem, env_name) {
  override <- Sys.getenv(env_name, unset = "")
  if (nzchar(override)) {
    if (!file.exists(override))
      stop(env_name, " points to a missing file: ", override, call. = FALSE)
    return(normalizePath(override, winslash = "/", mustWork = TRUE))
  }
  filenames <- c(paste0(stem, ".rda"), paste0(stem, ".RData"),
                 paste0(stem, ".csv"))
  dirs <- unique(c(
    study_b_data_dir,
    file.path(study_b_package_root, "01_DATA", "raw_rda"),
    file.path(study_b_package_root, "01_DATA"),
    file.path(study_b_package_root, "DATA")
  ))
  direct <- unlist(lapply(dirs, function(d) file.path(d, filenames)),
                   use.names = FALSE)
  direct <- direct[file.exists(direct)]
  if (length(direct))
    return(normalizePath(direct[1], winslash = "/", mustWork = TRUE))
  hits <- list.files(
    study_b_package_root,
    pattern = paste0("^", stem, "([^.]*)?[.](rda|RData|csv)$"),
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
  if (!length(hits)) return(NA_character_)
  hits <- hits[order(nchar(hits), hits)]
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}
.table1_read_data_file <- function(path, preferred_object) {
  if (grepl("[.]csv$", path, ignore.case = TRUE)) {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    sb_load_rda(path, preferred_object)
  }
}
.table1_pick_column <- function(dat, candidates, source_name, field_name) {
  hit <- candidates[candidates %in% names(dat)]
  if (!length(hit)) {
    stop(source_name, " has no recognized ", field_name, " column. Tried: ",
         paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit[1]
}
anchor_check <- sb_anchor %>%
  group_by(RID) %>%
  summarise(
    n_anchor_samples = n_distinct(sample_id),
    n_anchor_dates = n_distinct(plasma_date),
    n_pet_values = n_distinct(amyloid_centiloid, na.rm = TRUE),
    n_pet_classes = n_distinct(amyloid_positive, na.rm = TRUE),
    .groups = "drop"
  )
if (any(anchor_check$n_anchor_samples != 1L) ||
    any(anchor_check$n_anchor_dates != 1L) ||
    any(anchor_check$n_pet_values > 1L) ||
    any(anchor_check$n_pet_classes > 1L)) {
  stop("sb_anchor is not internally consistent within RID.", call. = FALSE)
}
if (!"phase" %in% names(sb_anchor)) {
  stop("sb_anchor has no `phase` (ORIGPROT) column; cannot report enrollment protocol.",
       call. = FALSE)
}
table1_phase_levels <- c("ADNI1", "ADNIGO", "ADNI2", "ADNI3", "ADNI4")
phase_check <- sb_anchor %>%
  group_by(RID) %>%
  summarise(n_phase = n_distinct(as.character(phase)), .groups = "drop")
if (any(phase_check$n_phase != 1L)) {
  stop("ORIGPROT (`phase`) is not constant within RID in sb_anchor.", call. = FALSE)
}
unrecognized_phase <- setdiff(unique(as.character(sb_anchor$phase)), table1_phase_levels)
if (length(unrecognized_phase)) {
  stop("Unrecognized ORIGPROT value(s): ", paste(unrecognized_phase, collapse = ", "),
       call. = FALSE)
}
anchor_person <- sb_anchor %>%
  arrange(RID, assay) %>%
  group_by(RID) %>%
  summarise(
    sample_id = first(sample_id),
    plasma_date = first(plasma_date),
    enroll_phase = first(as.character(phase)),
    amyloid_centiloid = sb_first_nonmissing(amyloid_centiloid),
    amyloid_positive = sb_first_nonmissing(amyloid_positive),
    amyloid_gap_days = sb_first_nonmissing(amyloid_gap_days),
    abs_amyloid_gap_days = sb_first_nonmissing(abs_amyloid_gap_days),
    .groups = "drop"
  ) %>%
  mutate(
    amyloid_group = factor(
      case_when(
        amyloid_positive == 0L ~ "Amyloid PET negative",
        amyloid_positive == 1L ~ "Amyloid PET positive",
        TRUE                   ~ "No eligible anchor PET"
      ),
      levels = c("Amyloid PET negative", "Amyloid PET positive",
                 "No eligible anchor PET")
    )
  )
.table1_expected <- c(
  participants = cohort$primary_anchor_participants,
  pet_available = cohort$anchor_events_with_pet,
  pet_negative = sum(sb_anchor$amyloid_positive == 0L, na.rm = TRUE) /
    length(study_b_assays),
  pet_positive = sum(sb_anchor$amyloid_positive == 1L, na.rm = TRUE) /
    length(study_b_assays)
)
.table1_observed <- c(
  participants = nrow(anchor_person),
  pet_available = sum(!is.na(anchor_person$amyloid_positive)),
  pet_negative = sum(anchor_person$amyloid_positive == 0L, na.rm = TRUE),
  pet_positive = sum(anchor_person$amyloid_positive == 1L, na.rm = TRUE)
)
if (!isTRUE(all.equal(unname(.table1_observed),
                      unname(.table1_expected), tolerance = 0))) {
  stop("Table 1 anchor counts do not match 02_core_results.rds. Expected ",
       paste(names(.table1_expected), .table1_expected, sep = "=", collapse = ", "),
       "; observed ",
       paste(names(.table1_observed), .table1_observed, sep = "=", collapse = ", "),
       ".", call. = FALSE)
}
demo_candidates <- anchor_person %>%
  select(RID, plasma_date) %>%
  left_join(
    ptdemog %>%
      transmute(
        RID = as.integer(as.character(RID)),
        demo_date = sb_safe_date(VISDATE),
        sex_raw = trimws(as.character(PTGENDER)),
        birth_month_year = trimws(as.character(PTDOB)),
        education = sb_safe_num(PTEDUCAT),
        race_raw = trimws(as.character(PTRACCAT)),
        ethnicity_raw = trimws(as.character(PTETHCAT))
      ),
    by = "RID", relationship = "many-to-many"
  ) %>%
  mutate(demo_gap_days = abs(as.numeric(demo_date - plasma_date))) %>%
  arrange(RID, is.na(demo_gap_days), demo_gap_days, demo_date)
demo_person <- demo_candidates %>%
  group_by(RID) %>%
  summarise(
    plasma_date = first(plasma_date),
    demographic_date = sb_first_nonmissing(demo_date),
    demographic_gap_days = sb_first_nonmissing(demo_gap_days),
    sex_raw = sb_first_nonmissing(na_if(sex_raw, "")),
    birth_month_year = sb_first_nonmissing(na_if(birth_month_year, "")),
    education = sb_first_nonmissing(education),
    race_raw = sb_first_nonmissing(na_if(race_raw, "")),
    ethnicity_raw = sb_first_nonmissing(na_if(ethnicity_raw, "")),
    .groups = "drop"
  ) %>%
  mutate(
    birth_date_approx = as.Date(paste0(birth_month_year, "/15"),
                                format = "%m/%Y/%d"),
    age_years = as.numeric(plasma_date - birth_date_approx) / 365.2425,
    female = sb_map_sex(sex_raw),
    race = sb_map_race(race_raw),
    hispanic_latino = sb_map_ethnicity(ethnicity_raw)
  )
apoe_long <- apoeres %>%
  transmute(RID = as.integer(as.character(RID)),
            apoe_genotype = gsub("|", "/", sb_code_label(GENOTYPE), fixed = TRUE)) %>%
  mutate(apoe_genotype = vapply(strsplit(apoe_genotype, "/", fixed = TRUE),
                                function(z) if (all(is.na(z))) NA_character_ else paste(sort(z), collapse = "/"), character(1))) %>%
  filter(RID %in% anchor_person$RID, !is.na(apoe_genotype)) %>%
  distinct(RID, apoe_genotype)
apoe_conflicts <- apoe_long %>% count(RID, name = "distinct_genotypes") %>%
  filter(distinct_genotypes > 1L)
if (nrow(apoe_conflicts)) {
  stop("Conflicting APOE genotypes for ", nrow(apoe_conflicts),
       " Table 1 participant(s).", call. = FALSE)
}
if (any(!grepl("^[234][/|][234]$", apoe_long$apoe_genotype)))
  stop("Unrecognized APOE genotype; check 2/3, 3/4 etc. in APOERES before reporting.")
apoe_person <- apoe_long %>%
  mutate(apoe4_carrier = grepl("4", apoe_genotype, fixed = TRUE))
dx_source <- dxsum %>% transmute(
  RID = as.integer(as.character(RID)), diagnosis_date = sb_safe_date(EXAMDATE),
  diagnosis_raw = trimws(as.character(DIAGNOSIS)),
  diagnosis = sb_map_diagnosis(diagnosis_raw)) %>%
  filter(RID %in% anchor_person$RID)
dx_conflicts <- dx_source %>%
  filter(!is.na(diagnosis), !is.na(diagnosis_date)) %>%
  group_by(RID, diagnosis_date) %>%
  summarise(
    n_diagnoses = n_distinct(diagnosis),
    raw_codes = paste(sort(unique(diagnosis_raw)), collapse = " | "),
    diagnoses = paste(sort(unique(diagnosis)), collapse = " | "),
    .groups = "drop"
  ) %>%
  filter(n_diagnoses > 1L)
dx_candidates_all <- anchor_person %>%
  select(RID, plasma_date) %>%
  inner_join(dx_source, by = "RID", relationship = "many-to-many") %>%
  mutate(
    diagnosis_gap_days = abs(as.numeric(diagnosis_date - plasma_date))
  ) %>%
  filter(
    is.finite(diagnosis_gap_days),
    diagnosis_gap_days <= 365,
    !is.na(diagnosis)
  ) %>%
  arrange(RID, diagnosis_gap_days, diagnosis_date)
dx_nearest_audit <- dx_candidates_all %>%
  group_by(RID) %>%
  filter(diagnosis_gap_days == min(diagnosis_gap_days)) %>%
  summarise(
    nearest_dates = n_distinct(diagnosis_date),
    nearest_diagnoses = n_distinct(diagnosis),
    diagnoses = paste(sort(unique(diagnosis)), collapse = " | "),
    raw_codes = paste(sort(unique(diagnosis_raw)), collapse = " | "),
    .groups = "drop"
  ) %>%
  filter(nearest_diagnoses > 1L)
if (nrow(dx_nearest_audit)) {
  print(dx_nearest_audit, n = Inf)
  stop(
    "Unresolved nearest-anchor DXSUM diagnosis conflict detected. ",
    "Review the printed RID(s) before generating Table 1.",
    call. = FALSE
  )
}
if (nrow(dx_conflicts)) {
  message(
    "DXSUM audit: ", nrow(dx_conflicts),
    " same-day diagnosis conflict cell(s) detected; none affect the ",
    "nearest-anchor diagnosis selection."
  )
}
unknown_dx <- dx_source %>% filter(is.na(diagnosis), !sb_declared_unknown("DIAGNOSIS", diagnosis_raw))
if (nrow(unknown_dx)) stop("Unrecognized DXSUM DIAGNOSIS code; inspect the source dictionary.")
dx_candidates <- anchor_person %>%
  select(RID, plasma_date) %>%
  inner_join(
    dx_source %>% select(-diagnosis),
    by = "RID", relationship = "many-to-many"
  ) %>%
  mutate(diagnosis_gap_days = abs(as.numeric(diagnosis_date - plasma_date))) %>%
  filter(is.finite(diagnosis_gap_days), diagnosis_gap_days <= 365,
         !is.na(sb_map_diagnosis(diagnosis_raw))) %>%
  arrange(RID, diagnosis_gap_days, diagnosis_date) %>%
  distinct(RID, .keep_all = TRUE) %>%
  mutate(
    diagnosis = sb_map_diagnosis(diagnosis_raw)
  ) %>%
  select(RID, diagnosis_date, diagnosis_gap_days, diagnosis_raw, diagnosis)
adnimerge_path <- .table1_find_data_file("ADNIMERGE", "STUDY_B_ADNIMERGE_FILE")
if (!is.na(adnimerge_path)) {
  adnimerge <- .table1_read_data_file(adnimerge_path, "ADNIMERGE")
  needed <- c("RID", "EXAMDATE", "MMSE", "CDRSB")
  missing <- setdiff(needed, names(adnimerge))
  if (length(missing))
    stop("ADNIMERGE is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  cognition_long <- adnimerge %>%
    transmute(
      RID = as.integer(as.character(RID)),
      cognition_date = sb_safe_date(EXAMDATE),
      MMSE = sb_safe_num(MMSE),
      CDRSB = sb_safe_num(CDRSB)
    )
} else {
  mmse_path <- .table1_find_data_file("MMSE", "STUDY_B_MMSE_FILE")
  cdr_path <- .table1_find_data_file("CDR", "STUDY_B_CDR_FILE")
  if (is.na(mmse_path) || is.na(cdr_path)) {
    stop(
      "Cannot locate cognitive data. Supply ADNIMERGE.rda/.csv, or both ",
      "MMSE.rda/.csv and CDR.rda/.csv under STUDY_B_ROOT. You may also set ",
      "STUDY_B_ADNIMERGE_FILE, STUDY_B_MMSE_FILE and STUDY_B_CDR_FILE.",
      call. = FALSE
    )
  }
  mmse_raw <- .table1_read_data_file(mmse_path, "MMSE")
  cdr_raw <- .table1_read_data_file(cdr_path, "CDR")
  if (!"RID" %in% names(mmse_raw) || !"RID" %in% names(cdr_raw))
    stop("MMSE and CDR files must both contain RID.", call. = FALSE)
  mmse_date_col <- .table1_pick_column(
    mmse_raw, c("EXAMDATE", "VISDATE"), basename(mmse_path), "date")
  mmse_score_col <- .table1_pick_column(
    mmse_raw, c("MMSCORE", "MMSE"), basename(mmse_path), "MMSE score")
  mmse_long <- tibble(
    RID = as.integer(as.character(mmse_raw[["RID"]])),
    cognition_date = sb_safe_date(mmse_raw[[mmse_date_col]]),
    MMSE = sb_safe_num(mmse_raw[[mmse_score_col]])
  )
  cdr_date_col <- .table1_pick_column(
    cdr_raw, c("EXAMDATE", "VISDATE"), basename(cdr_path), "date")
  if ("CDRSB" %in% names(cdr_raw)) {
    cdrsb_value <- sb_safe_num(cdr_raw[["CDRSB"]])
  } else {
    cdr_boxes <- c("CDMEMORY", "CDORIENT", "CDJUDGE",
                   "CDCOMMUN", "CDHOME", "CDCARE")
    missing_boxes <- setdiff(cdr_boxes, names(cdr_raw))
    if (length(missing_boxes)) {
      stop("CDR has neither CDRSB nor all six box fields; missing: ",
           paste(missing_boxes, collapse = ", "), call. = FALSE)
    }
    box_matrix <- as.data.frame(lapply(cdr_raw[cdr_boxes], sb_safe_num))
    cdrsb_value <- ifelse(rowSums(is.finite(as.matrix(box_matrix))) == 6L,
                          rowSums(box_matrix), NA_real_)
  }
  cdr_long <- tibble(
    RID = as.integer(as.character(cdr_raw[["RID"]])),
    cognition_date = sb_safe_date(cdr_raw[[cdr_date_col]]),
    CDRSB = cdrsb_value
  )
  cognition_long <- full_join(mmse_long, cdr_long,
                              by = c("RID", "cognition_date"),
                              relationship = "many-to-many")
}
cognition_long <- cognition_long %>%
  filter(RID %in% anchor_person$RID, !is.na(cognition_date))
cognition_conflicts <- cognition_long %>%
  group_by(RID, cognition_date) %>%
  summarise(
    n_mmse = n_distinct(MMSE[is.finite(MMSE)]),
    n_cdrsb = n_distinct(CDRSB[is.finite(CDRSB)]),
    .groups = "drop"
  ) %>%
  filter(n_mmse > 1L | n_cdrsb > 1L)
if (nrow(cognition_conflicts)) {
  stop("Conflicting MMSE or CDR-SB values in ", nrow(cognition_conflicts),
       " RID-date cell(s) in the cognitive source data.", call. = FALSE)
}
mmse_person <- anchor_person %>%
  select(RID, plasma_date) %>%
  inner_join(cognition_long %>% filter(is.finite(MMSE), between(MMSE, 0, 30)),
             by = "RID", relationship = "many-to-many") %>%
  mutate(mmse_gap_days = abs(as.numeric(cognition_date - plasma_date))) %>%
  filter(mmse_gap_days <= 365) %>%
  arrange(RID, mmse_gap_days, cognition_date) %>%
  distinct(RID, .keep_all = TRUE) %>%
  transmute(RID, mmse = MMSE, mmse_date = cognition_date, mmse_gap_days)
cdrsb_person <- anchor_person %>%
  select(RID, plasma_date) %>%
  inner_join(cognition_long %>% filter(is.finite(CDRSB), between(CDRSB, 0, 18)),
             by = "RID", relationship = "many-to-many") %>%
  mutate(cdrsb_gap_days = abs(as.numeric(cognition_date - plasma_date))) %>%
  filter(cdrsb_gap_days <= 365) %>%
  arrange(RID, cdrsb_gap_days, cognition_date) %>%
  distinct(RID, .keep_all = TRUE) %>%
  transmute(RID, cdrsb = CDRSB, cdrsb_date = cognition_date, cdrsb_gap_days)
followup_person <- sb_events %>%
  distinct(RID, sample_id, plasma_date) %>%
  arrange(RID, plasma_date) %>%
  group_by(RID) %>%
  summarise(
    collection_occasions = n(),
    followup_years = as.numeric(last(plasma_date) - first(plasma_date)) / 365.2425,
    .groups = "drop"
  )
assay_person <- sb_anchor %>%
  transmute(RID, assay = as.character(assay), raw_value = exp(log_value)) %>%
  tidyr::pivot_wider(names_from = assay, values_from = raw_value,
                     names_prefix = "assay_")
table1_person <- anchor_person %>%
  left_join(select(demo_person, -plasma_date), by = "RID") %>%
  left_join(apoe_person, by = "RID") %>%
  left_join(dx_candidates, by = "RID") %>%
  left_join(mmse_person, by = "RID") %>%
  left_join(cdrsb_person, by = "RID") %>%
  left_join(followup_person, by = "RID") %>%
  left_join(assay_person, by = "RID")
if (nrow(table1_person) != n_distinct(sb_anchor$RID) ||
    anyDuplicated(table1_person$RID)) {
  stop("Table 1 participant dataset is not one row per anchor RID.", call. = FALSE)
}
fmt_mean_sd <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("--")
  sprintf("%.1f (%.1f)", mean(x), stats::sd(x))
}
fmt_median_iqr <- function(x, digits = 1L) {
  x <- x[is.finite(x)]
  if (!length(x)) return("--")
  q <- stats::quantile(x, c(.25, .50, .75), na.rm = TRUE, names = FALSE)
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
          q[2], q[1], q[3])
}
fmt_n_pct <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("--")
  sprintf("%d (%.1f%%)", sum(x), 100 * mean(x))
}
table1_specs <- list(
  list("Age at anchor, years", "Mean (SD)",
       function(d) fmt_mean_sd(d$age_years), function(d) sum(is.finite(d$age_years))),
  list("Female sex", "n (%)", function(d) fmt_n_pct(d$female),
       function(d) sum(!is.na(d$female))),
  list("Education, years", "Mean (SD)",
       function(d) fmt_mean_sd(d$education), function(d) sum(is.finite(d$education))),
  list("APOE e4 carrier", "n (%)", function(d) fmt_n_pct(d$apoe4_carrier),
       function(d) sum(!is.na(d$apoe4_carrier))),
  list("White race", "n (%)", function(d) fmt_n_pct(d$race == "White"),
       function(d) sum(!is.na(d$race))),
  list("Black or African American race", "n (%)",
       function(d) fmt_n_pct(d$race == "Black or African American"),
       function(d) sum(!is.na(d$race))),
  list("Asian race", "n (%)", function(d) fmt_n_pct(d$race == "Asian"),
       function(d) sum(!is.na(d$race))),
  list("Other or multiple race", "n (%)",
       function(d) fmt_n_pct(d$race == "Other or multiple race"),
       function(d) sum(!is.na(d$race))),
  list("Hispanic or Latino ethnicity", "n (%)",
       function(d) fmt_n_pct(d$hispanic_latino),
       function(d) sum(!is.na(d$hispanic_latino))),
  list("Cognitively unimpaired", "n (%)",
       function(d) fmt_n_pct(d$diagnosis == "Cognitively unimpaired"),
       function(d) sum(!is.na(d$diagnosis))),
  list("Mild cognitive impairment", "n (%)",
       function(d) fmt_n_pct(d$diagnosis == "MCI"),
       function(d) sum(!is.na(d$diagnosis))),
  list("Dementia", "n (%)", function(d) fmt_n_pct(d$diagnosis == "Dementia"),
       function(d) sum(!is.na(d$diagnosis))),
  list("ADNI enrollment protocol: ADNI1", "n (%)",
       function(d) fmt_n_pct(d$enroll_phase == "ADNI1"), function(d) sum(!is.na(d$enroll_phase))),
  list("ADNI enrollment protocol: ADNIGO", "n (%)",
       function(d) fmt_n_pct(d$enroll_phase == "ADNIGO"), function(d) sum(!is.na(d$enroll_phase))),
  list("ADNI enrollment protocol: ADNI2", "n (%)",
       function(d) fmt_n_pct(d$enroll_phase == "ADNI2"), function(d) sum(!is.na(d$enroll_phase))),
  list("ADNI enrollment protocol: ADNI3", "n (%)",
       function(d) fmt_n_pct(d$enroll_phase == "ADNI3"), function(d) sum(!is.na(d$enroll_phase))),
  list("ADNI enrollment protocol: ADNI4", "n (%)",
       function(d) fmt_n_pct(d$enroll_phase == "ADNI4"), function(d) sum(!is.na(d$enroll_phase))),
  list("MMSE score", "Mean (SD)",
       function(d) fmt_mean_sd(d$mmse), function(d) sum(is.finite(d$mmse))),
  list("CDR-SB score", "Median [IQR]",
       function(d) fmt_median_iqr(d$cdrsb, 1L),
       function(d) sum(is.finite(d$cdrsb))),
  list("Amyloid PET reference available at anchor", "n (%)",
       function(d) fmt_n_pct(!is.na(d$amyloid_positive)), function(d) nrow(d)),
  list("Centiloid at anchor", "Median [IQR]",
       function(d) fmt_median_iqr(d$amyloid_centiloid, 1L),
       function(d) sum(is.finite(d$amyloid_centiloid))),
  list("Absolute plasma-PET interval, days", "Median [IQR]",
       function(d) fmt_median_iqr(d$abs_amyloid_gap_days, 0L),
       function(d) sum(is.finite(d$abs_amyloid_gap_days))),
  list("Collection occasions per participant", "Median [IQR]",
       function(d) fmt_median_iqr(d$collection_occasions, 0L),
       function(d) sum(is.finite(d$collection_occasions))),
  list("First-to-last collection interval, years", "Median [IQR]",
       function(d) fmt_median_iqr(d$followup_years, 2L),
       function(d) sum(is.finite(d$followup_years))),
  list("ALZpath p-tau217, pg/mL", "Median [IQR]",
       function(d) fmt_median_iqr(d$assay_AlzPath, 3L),
       function(d) sum(is.finite(d$assay_AlzPath))),
  list("C2N %p-tau217, %", "Median [IQR]",
       function(d) fmt_median_iqr(d$assay_C2N_ratio, 3L),
       function(d) sum(is.finite(d$assay_C2N_ratio))),
  list("Fujirebio p-tau217, pg/mL", "Median [IQR]",
       function(d) fmt_median_iqr(d$assay_Fujirebio, 3L),
       function(d) sum(is.finite(d$assay_Fujirebio))),
  list("Janssen p-tau217, pg/mL", "Median [IQR]",
       function(d) fmt_median_iqr(d$assay_Janssen, 3L),
       function(d) sum(is.finite(d$assay_Janssen)))
)
group_levels <- c("Overall", levels(table1_person$amyloid_group))
group_data <- c(
  list(Overall = table1_person),
  setNames(lapply(levels(table1_person$amyloid_group), function(g) {
    table1_person[table1_person$amyloid_group == g, , drop = FALSE]
  }), levels(table1_person$amyloid_group))
)
group_labels <- setNames(
  sprintf("%s (N=%d)", group_levels,
          vapply(group_data[group_levels], nrow, integer(1))),
  group_levels
)
.present_phases <- unique(table1_person$enroll_phase)
table1_specs <- Filter(function(s) {
  !startsWith(s[[1]], "ADNI enrollment protocol: ") ||
    sub("^ADNI enrollment protocol: ", "", s[[1]]) %in% .present_phases
}, table1_specs)
table1_long <- purrr::map_dfr(seq_along(table1_specs), function(i) {
  s <- table1_specs[[i]]
  purrr::map_dfr(group_levels, function(g) {
    d <- group_data[[g]]
    tibble(
      row_order = i,
      characteristic = s[[1]],
      statistic = s[[2]],
      group = unname(group_labels[[g]]),
      value = s[[3]](d),
      n_nonmissing = as.integer(s[[4]](d)),
      group_n = nrow(d)
    )
  })
})
table1_formatted <- table1_long %>%
  select(row_order, characteristic, statistic, group, value) %>%
  pivot_wider(names_from = group, values_from = value) %>%
  arrange(row_order) %>%
  select(-row_order)
table1_denominators <- table1_long %>%
  select(characteristic, statistic, group, n_nonmissing, group_n)
table1_source_audit <- tibble(
  item = c(
    "Anchor participants",
    "Anchor participants with eligible amyloid PET",
    "Anchor amyloid PET negative",
    "Anchor amyloid PET positive",
    "Anchor participants without eligible amyloid PET",
    "Missing age",
    "Missing sex",
    "Missing education",
    "Missing APOE genotype",
    "Missing diagnosis within 365 days",
    "Missing MMSE within 365 days",
    "Missing CDR-SB within 365 days",
    "Conflicting APOE genotypes"
  ),
  n = c(
    nrow(table1_person),
    sum(!is.na(table1_person$amyloid_positive)),
    sum(table1_person$amyloid_positive == 0L, na.rm = TRUE),
    sum(table1_person$amyloid_positive == 1L, na.rm = TRUE),
    sum(is.na(table1_person$amyloid_positive)),
    sum(!is.finite(table1_person$age_years)),
    sum(is.na(table1_person$female)),
    sum(!is.finite(table1_person$education)),
    sum(is.na(table1_person$apoe_genotype)),
    sum(is.na(table1_person$diagnosis)),
    sum(!is.finite(table1_person$mmse)),
    sum(!is.finite(table1_person$cdrsb)),
    nrow(apoe_conflicts)
  )
)
table1_encoding_audit <- bind_rows(
  demo_person %>% transmute(field="PTGENDER",raw=sex_raw,derived=as.character(female)),
  demo_person %>% transmute(field="PTRACCAT",raw=race_raw,derived=race),
  demo_person %>% transmute(field="PTETHCAT",raw=ethnicity_raw,derived=as.character(hispanic_latino)),
  dx_candidates %>% transmute(field="DIAGNOSIS",raw=diagnosis_raw,derived=diagnosis)
) %>% count(field,raw,derived,name="participants")
sb_write(table1_encoding_audit,"Table_1_encoding_audit.csv")
unrecognized_encoding <- table1_encoding_audit %>%
  filter(is.na(derived), !sb_declared_unknown(field, raw))
if (nrow(unrecognized_encoding)) {
  print(unrecognized_encoding, n = Inf)
  stop("Unrecognized Table 1 field codes are printed above and saved in ",
       "Table_1_encoding_audit.csv. Do not report Table 1 until mapped.",
       call. = FALSE)
}
sb_write(table1_formatted,    "Table_1_participant_characteristics.csv")
sb_write(table1_denominators, "Table_1_denominators.csv")
sb_write(table1_person,       "Table_1_participant_audit.csv")
sb_write(table1_source_audit, "Table_1_source_audit.csv")
cat("\nTable 1 -- participant demographic and clinical characteristics\n")
print(as.data.frame(table1_formatted), row.names = FALSE)
cat("\nTable 1 -- source and missingness audit\n")
print(as.data.frame(table1_source_audit), row.names = FALSE)
.table1_used_inputs <- c(unname(table1_files),
                         if (!is.na(adnimerge_path)) adnimerge_path else c(mmse_path, cdr_path))
for (.file in c("Table_1_participant_characteristics.csv", "Table_1_denominators.csv",
                "Table_1_participant_audit.csv", "Table_1_source_audit.csv"))
  sb_stamp_result(file.path(study_b_dir, .file), .table1_used_inputs)
