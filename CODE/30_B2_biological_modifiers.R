options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(lme4)
})
ROOT <- if (exists("study_b_package_root")) study_b_package_root else Sys.getenv("STUDY_B_ROOT")
RES <- if (exists("study_b_dir")) study_b_dir else file.path(ROOT, "RESULTS")
RAW <- if (exists("study_b_data_dir")) study_b_data_dir else
  Sys.getenv("STUDY_B_DATA_DIR", unset = file.path(ROOT, "01_DATA", "raw_rda"))
ASSAYS <- c("AlzPath", "C2N_ratio", "Fujirebio", "Janssen")
B2_BOOT <- as.integer(Sys.getenv("STUDY_B_B2_BOOT", unset = "1000"))
B2_CORES <- as.integer(Sys.getenv("STUDY_B_B2_CORES", unset = "1"))
SEED <- 20260917L
stopifnot(!is.na(B2_BOOT), B2_BOOT >= 0L, !is.na(B2_CORES), B2_CORES >= 1L)
core <- sb_read_core(file.path(RES, "02_core_results.rds"))
persons <- sb_read_result(file.path(RES, "Table_1_participant_audit.csv"), check.names = FALSE)
persons <- persons %>% transmute(RID = as.integer(as.character(RID)),
  birth_date = as.Date(birth_date_approx), female = as.logical(female))
LABDATA <- sb_load_rda(file.path(RAW, "LABDATA.rda"), "LABDATA")
VITALS <- sb_load_rda(file.path(RAW, "VITALS.rda"), "VITALS")
for (nm in c("RID", "EXAMDATE", "RCT392", "HMT40", "HMT2", "HMT3", "RCT6"))
  if (!nm %in% names(LABDATA)) stop("LABDATA missing required column: ", nm)
for (nm in c("RID", "VSWEIGHT", "VSWTUNIT", "VSHEIGHT", "VSHTUNIT"))
  if (!nm %in% names(VITALS)) stop("VITALS missing required column: ", nm)
VITAL_DATE_FIELD <- if ("EXAMDATE" %in% names(VITALS)) "EXAMDATE" else
  if ("VISDATE" %in% names(VITALS)) "VISDATE" else
    stop("VITALS requires EXAMDATE (official field) or the documented legacy VISDATE.")
WEIGHT_CODE_MAP <- character()
HEIGHT_CODE_MAP <- character()
normalise_vital_unit <- function(x, map, kind) {
  u <- tolower(trimws(as.character(x)))
  u[u %in% c("", "na", "nan", "-1", "-4", "-9")] <- NA_character_
  if (length(map)) {
    if (is.null(names(map)) || anyDuplicated(names(map)))
      stop("Unit code maps must have unique code names.")
    hit <- !is.na(u) & u %in% names(map)
    u[hit] <- unname(map[u[hit]])
  }
  if (kind == "weight") {
    out <- case_when(u %in% c("lb", "lbs", "pound", "pounds") ~ "lb",
                     u %in% c("kg", "kilogram", "kilograms") ~ "kg",
                     TRUE ~ NA_character_)
  } else {
    out <- case_when(u %in% c("in", "inch", "inches") ~ "in",
                     u %in% c("cm", "centimeter", "centimeters", "centimetre", "centimetres") ~ "cm",
                     u %in% c("m", "meter", "meters", "metre", "metres") ~ "m",
                     TRUE ~ NA_character_)
  }
  bad <- unique(u[!is.na(u) & is.na(out)])
  if (length(bad)) stop("Unresolved VITALS ", kind, " unit(s): ",
    paste(bad, collapse = ", "), ". Set ", toupper(kind),
    "_CODE_MAP from your release dictionary; do not guess units.")
  out
}
if (anyDuplicated(persons$RID)) stop("B2 demographics must have one row per RID.")
num <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))
date_safe <- sb_safe_date
CREAT_UNIT <- tolower(Sys.getenv("STUDY_B_CREATININE_UNIT", unset = "mg/dL"))
HB_UNIT <- tolower(Sys.getenv("STUDY_B_HEMOGLOBIN_UNIT", unset = "g/dL"))
stopifnot(CREAT_UNIT %in% c("mg/dl", "umol/l", "µmol/l"),
          HB_UNIT %in% c("g/dl", "g/l"))
creat_num <- function(x) num(x) / if (CREAT_UNIT == "mg/dl") 1 else 88.4
hb_num <- function(x) num(x) / if (HB_UNIT == "g/dl") 1 else 10
labs <- LABDATA %>% transmute(
  RID = as.integer(as.character(RID)), date = date_safe(EXAMDATE),
  creatinine_mg_dl = creat_num(RCT392), hemoglobin_g_dl = hb_num(HMT40),
  hematocrit_pct = num(HMT2), rbc_million_ul = num(HMT3), bun_mg_dl = num(RCT6)) %>%
  filter(!is.na(RID), !is.na(date)) %>%
  mutate(creatinine_valid = between(creatinine_mg_dl, 0.2, 15),
         hemoglobin_valid = between(hemoglobin_g_dl, 3, 25),
         hematocrit_valid = between(hematocrit_pct, 10, 75),
         rbc_valid = between(rbc_million_ul, 1, 10),
         bun_valid = between(bun_mg_dl, 1, 150),
         creatinine_mg_dl = if_else(creatinine_valid, creatinine_mg_dl, NA_real_),
         hemoglobin_g_dl = if_else(hemoglobin_valid, hemoglobin_g_dl, NA_real_),
         hematocrit_pct = if_else(hematocrit_valid, hematocrit_pct, NA_real_),
         rbc_million_ul = if_else(rbc_valid, rbc_million_ul, NA_real_),
         bun_mg_dl = if_else(bun_valid, bun_mg_dl, NA_real_))
vitals <- VITALS %>% transmute(
  RID = as.integer(as.character(RID)), date = date_safe(.data[[VITAL_DATE_FIELD]]),
  weight = num(VSWEIGHT),
  weight_unit = normalise_vital_unit(VSWTUNIT, WEIGHT_CODE_MAP, "weight"),
  height = num(VSHEIGHT),
  height_unit = normalise_vital_unit(VSHTUNIT, HEIGHT_CODE_MAP, "height")) %>%
  mutate(weight_kg = case_when(weight_unit == "lb" ~ weight*0.45359237,
                               weight_unit == "kg" ~ weight,
                               TRUE ~ NA_real_),
         height_m = case_when(height_unit == "in" ~ height*0.0254,
                              height_unit == "cm" ~ height/100,
                              height_unit == "m" ~ height,
                              TRUE ~ NA_real_),
         measurement_valid = between(weight_kg, 25, 300) & between(height_m, 1.2, 2.3),
         bmi = weight_kg/height_m^2,
         bmi_valid = measurement_valid & between(bmi, 10, 80),
         bmi = if_else(bmi_valid, bmi, NA_real_)) %>%
  filter(!is.na(RID), !is.na(date))
sb_save_csv(vitals, file.path(RES, "Table_S11_B2_vitals_unit_audit.csv"),
            row.names = FALSE, na = "")
events <- core$sb_events %>% transmute(RID = as.integer(as.character(RID)),
  sample_id = as.integer(as.character(sample_id)), plasma_date = as.Date(plasma_date),
  primary_anchor = as.logical(primary_anchor))
stopifnot("B2 requires globally unique sample_id in sb_events" =
            !anyDuplicated(events$sample_id))
nearest <- function(events, source, value_cols, prefix) {
  source <- source %>% filter(RID %in% events$RID)
  conflict <- source %>% group_by(RID, date) %>%
    summarise(across(all_of(value_cols), ~n_distinct(.x[is.finite(.x)])), .groups = "drop") %>%
    filter(if_any(all_of(value_cols), ~.x > 1L))
  if (nrow(conflict)) stop("Conflicting ", prefix,
    "values for the same RID/date; resolve source records before nearest-date matching.")
  events %>%
    inner_join(source %>% select(RID, source_date = date, all_of(value_cols)),
               by = "RID", relationship = "many-to-many") %>%
    mutate(gap_days = abs(as.numeric(source_date-plasma_date))) %>%
    arrange(sample_id, gap_days, source_date) %>% distinct(sample_id, .keep_all = TRUE) %>%
    select(sample_id, source_date, gap_days, all_of(value_cols)) %>%
    rename_with(~paste0(prefix, .x), -sample_id)
}
creat_match <- nearest(events, filter(labs, is.finite(creatinine_mg_dl)),
                       c("creatinine_mg_dl"), "creat_")
hb_match <- nearest(events, filter(labs, is.finite(hemoglobin_g_dl)),
                    c("hemoglobin_g_dl"), "hb_")
vital_match <- nearest(events, filter(vitals, is.finite(bmi)), c("bmi"), "vital_")
matched <- events %>% left_join(creat_match, by = "sample_id") %>%
  left_join(hb_match, by = "sample_id") %>%
  left_join(vital_match, by = "sample_id") %>% left_join(persons, by = "RID") %>%
  mutate(age_at_lab = as.numeric(creat_source_date-birth_date)/365.2425,
         kappa = if_else(female, 0.7, 0.9),
         alpha = if_else(female, -0.241, -0.302),
         scr_k = creat_creatinine_mg_dl/kappa,
         egfr = 142*pmin(scr_k, 1)^alpha*pmax(scr_k, 1)^(-1.200)*
           0.9938^age_at_lab*if_else(female, 1.012, 1),
         egfr = if_else(is.finite(egfr) & between(age_at_lab, 18, 120) &
                          between(egfr, 1, 200), egfr, NA_real_))
sb_save_csv(matched, file.path(RES, "Table_S11_B2_matched_measurement_audit.csv"),
            row.names = FALSE, na = "")
sb_save_csv(tibble(field = c("RCT392", "HMT40"),
                   input_unit = c(CREAT_UNIT, HB_UNIT),
                   analysis_unit = c("mg/dL", "g/dL"),
                   note = "Release-specific unit assumption; verify against the source laboratory documentation"),
            file.path(RES, "Table_S11_B2_unit_contract.csv"), row.names = FALSE)
coverage_one <- function(d, scope, modifier, value, gap) {
  value <- rlang::ensym(value); gap <- rlang::ensym(gap)
  map_dfr(c(0, 30, 90, 180, 365), function(w) tibble(
    scope = scope, modifier = modifier, window_days = w,
    denominator = nrow(d),
    available_n = sum(is.finite(dplyr::pull(d, !!value)) &
                        dplyr::pull(d, !!gap) <= w, na.rm = TRUE),
    available_pct = 100*available_n/denominator,
    median_abs_gap_days = if (available_n) median(dplyr::pull(d, !!gap)[
      is.finite(dplyr::pull(d, !!value)) & dplyr::pull(d, !!gap) <= w], na.rm=TRUE) else NA_real_))
}
coverage <- imap_dfr(list("Primary anchor"=filter(matched, primary_anchor %in% TRUE),
                          "All occasions"=matched), function(d, scope) {
  bind_rows(
    coverage_one(d, scope, "Creatinine", creat_creatinine_mg_dl, creat_gap_days),
    coverage_one(d, scope, "eGFR", egfr, creat_gap_days),
    coverage_one(d, scope, "Hemoglobin", hb_hemoglobin_g_dl, hb_gap_days),
    coverage_one(d, scope, "BMI", vital_bmi, vital_gap_days))
})
sb_save_csv(coverage, file.path(RES, "Table_S11_B2_modifier_coverage.csv"), row.names=FALSE)
value_audit <- bind_rows(
  tibble(modifier="Creatinine, mg/dL", source="LABDATA RCT392",
         source_nonmissing=sum(is.finite(creat_num(LABDATA$RCT392))),
         source_implausible=sum(is.finite(creat_num(LABDATA$RCT392)) & !between(creat_num(LABDATA$RCT392),.2,15)),
         minimum=min(labs$creatinine_mg_dl,na.rm=TRUE), median=median(labs$creatinine_mg_dl,na.rm=TRUE),
         maximum=max(labs$creatinine_mg_dl,na.rm=TRUE)),
  tibble(modifier="Hemoglobin, g/dL", source="LABDATA HMT40",
         source_nonmissing=sum(is.finite(hb_num(LABDATA$HMT40))),
         source_implausible=sum(is.finite(hb_num(LABDATA$HMT40)) & !between(hb_num(LABDATA$HMT40),3,25)),
         minimum=min(labs$hemoglobin_g_dl,na.rm=TRUE), median=median(labs$hemoglobin_g_dl,na.rm=TRUE),
         maximum=max(labs$hemoglobin_g_dl,na.rm=TRUE)),
  tibble(modifier="BMI, kg/m2", source="VITALS explicit fields/units",
         source_nonmissing=sum(is.finite(vitals$bmi)),
         source_implausible=sum(!vitals$bmi_valid & is.finite(vitals$weight_kg) &
                                  is.finite(vitals$height_m), na.rm=TRUE),
         minimum=min(vitals$bmi,na.rm=TRUE), median=median(vitals$bmi,na.rm=TRUE), maximum=max(vitals$bmi,na.rm=TRUE)))
sb_save_csv(value_audit, file.path(RES, "Table_S11_B2_modifier_value_audit.csv"), row.names=FALSE)
mods <- matched %>% filter(primary_anchor) %>%
  mutate(egfr_90 = if_else(creat_gap_days <= 90, egfr, NA_real_),
         hb_90 = if_else(hb_gap_days <= 90, hb_hemoglobin_g_dl, NA_real_),
         bmi_180 = if_else(vital_gap_days <= 180, vital_bmi, NA_real_)) %>%
  select(RID, egfr_90, hb_90, bmi_180)
stopifnot(!anyDuplicated(mods$RID))
mods <- mods %>% filter(if_all(c(egfr_90, hb_90, bmi_180), is.finite))
for (nm in c("egfr_90", "hb_90", "bmi_180")) {
  z <- mods[[nm]]
  if (length(z) > 1L && is.finite(sd(z)) && sd(z) > 0)
    mods[[paste0(nm, "_z")]] <- as.numeric(scale(z))
  else mods[[paste0(nm, "_z")]] <- NA_real_
}
d <- core$sb_decision_long %>%
  mutate(RID = as.integer(as.character(RID)), sample_id=as.integer(as.character(sample_id))) %>%
  inner_join(mods, by="RID") %>%
  filter(is.finite(egfr_90), is.finite(hb_90), is.finite(bmi_180)) %>%
  mutate(RID=factor(RID), assay=factor(assay,levels=ASSAYS), phase=factor(phase),
         occ_id=interaction(RID,sample_id,drop=TRUE),
         baseline_amyloid_group=factor(baseline_amyloid_group))
MODEL_SPECS <- list(
  same_subset_main_effects_only = character(),
  plus_assay_by_eGFR            = "assay:egfr_90_z",
  plus_assay_by_hemoglobin      = "assay:hb_90_z",
  plus_assay_by_BMI             = "assay:bmi_180_z",
  plus_all_assay_by_modifiers   = c("assay:egfr_90_z","assay:hb_90_z","assay:bmi_180_z"))
model_formula <- function(interactions=character(), data=d) {
  base <- c(sb_vc_fixed_terms(data),
            "egfr_90_z", "hb_90_z", "bmi_180_z")
  rhs <- paste(c(base, interactions, "(1|RID)", "(1|occ_id)", "(1|RID:assay)"), collapse=" + ")
  as.formula(paste("margin_youden ~", rhs))
}
fit_model <- function(label, interactions=character()) {
  m <- sb_checked_lmer(model_formula(interactions, data=d), data=d, REML=TRUE,
            control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  vc <- as.data.frame(VarCorr(m)); gv <- function(g) vc$vcov[match(g,vc$grp)]
  tibble(model=label,n_participants=n_distinct(d$RID),n_occasions=n_distinct(d$occ_id),n_rows=nrow(d),
         person_by_assay_variance=gv("RID:assay"),shared_occasion_variance=gv("occ_id"),
         residual_variance=gv("Residual"),ratio_person_assay_to_occasion=sb_safe_ratio(person_by_assay_variance, shared_occasion_variance),
         singular=isSingular(m),optimizer_code=if(is.null(m@optinfo$conv$opt))0L else m@optinfo$conv$opt[1],
         convergence_message=paste(m@optinfo$conv$lme4$messages,collapse="; "))
}
repeated_ids <- d %>% distinct(RID, sample_id) %>% count(RID) %>% filter(n >= 2L) %>% pull(RID)
model_identifiable <- n_distinct(d$RID) >= 100L && length(repeated_ids) >= 2L &&
  all(vapply(d[c("egfr_90_z", "hb_90_z", "bmi_180_z")],
             function(x) all(is.finite(x)), logical(1)))
if (model_identifiable) {
  models <- imap_dfr(MODEL_SPECS, function(ix, lab) fit_model(lab, ix)) %>%
    mutate(reference_variance=first(person_by_assay_variance),
           fraction_remaining=if_else(reference_variance > 1e-10,
                                     person_by_assay_variance/reference_variance, NA_real_),
           percent_explained=100*(1-fraction_remaining))
  boot_rows <- split(seq_len(nrow(d)), as.character(d$RID))
  set.seed(SEED + 5000L)
  boot_ids <- replicate(B2_BOOT, sample(names(boot_rows), length(boot_rows), replace=TRUE),
                        simplify=FALSE)
  d_boot_src <- d %>% mutate(RID = as.character(RID), sample_id = as.character(sample_id))
  boot_one <- function(ids) {
    idx <- boot_rows[ids]
    bd <- as.data.frame(d_boot_src)[unlist(idx, use.names=FALSE), , drop=FALSE]
    bd$RID <- factor(rep(paste0("b", seq_along(ids)), lengths(idx)))
    bd$occ_id <- interaction(bd$RID, bd$sample_id, drop=TRUE)
    errors <- rep(NA_character_, length(MODEL_SPECS))
    singular <- rep(NA, length(MODEL_SPECS))
    k <- 0L
    out <- vapply(MODEL_SPECS, function(ix) {
      k <<- k + 1L
      v <- tryCatch({
        m <- suppressMessages(suppressWarnings(
          sb_checked_lmer(model_formula(ix, data=bd), data=bd, REML=TRUE,
                     control=lme4::lmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))))
        singular[k] <<- lme4::isSingular(m)
        vc <- as.data.frame(lme4::VarCorr(m))
        vc$vcov[match("RID:assay", vc$grp)]
      }, error=function(e) { errors[k] <<- conditionMessage(e); NA_real_ })
      v
    }, numeric(1))
    attr(out, "errors") <- errors
    attr(out, "singular") <- singular
    out
  }
  if (B2_BOOT > 0) {
    message("B2 bootstrap: ", B2_BOOT, " replicates x ", length(MODEL_SPECS),
            " models on ", B2_CORES, " core(s)")
    if (B2_CORES > 1) {
      cl <- parallel::makeCluster(B2_CORES)
      boot_var <- tryCatch({
        parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(lme4)))
        parallel::clusterExport(cl, c("d_boot_src", "boot_rows", "MODEL_SPECS", "model_formula",
                                  "sb_vc_fixed_terms", "sb_checked_lmer", "sb_lmer_converged",
                                  "sb_model_audit"),
                                envir=globalenv())
        parallel::parLapplyLB(cl, boot_ids, boot_one)
      }, finally = parallel::stopCluster(cl))
    } else {
      boot_var <- lapply(boot_ids, boot_one)
    }
    boot_diagnostics <- map_dfr(seq_along(boot_var), function(b) tibble(
      replicate = b, model = names(MODEL_SPECS),
      finite_variance = is.finite(boot_var[[b]]),
      singular = attr(boot_var[[b]], "singular"),
      error = attr(boot_var[[b]], "errors")))
    sb_save_csv(boot_diagnostics, file.path(RES, "Table_S11_B2_bootstrap_diagnostics.csv"),
                row.names = FALSE, na = "")
    boot_var <- do.call(rbind, boot_var)
    boot_ref <- boot_var[, 1]
    boot_pct <- 100 * (1 - boot_var / boot_ref)
    usable <- is.finite(boot_ref) & boot_ref > 1e-10
    boot_pct[!usable, ] <- NA_real_
    sb_save_csv(data.frame(replicate=seq_len(nrow(boot_var)), boot_var, check.names=FALSE),
              file.path(RES,"Table_S11_B2_bootstrap_replicates.csv"), row.names=FALSE, na="")
    boot_ci <- tibble(
      model = names(MODEL_SPECS),
      percent_explained_low  = apply(boot_pct, 2, function(z) sb_ci(z)[["low"]]),
      percent_explained_high = apply(boot_pct, 2, function(z) sb_ci(z)[["high"]]),
      person_by_assay_variance_low  = apply(boot_var, 2, function(z) sb_ci(z)[["low"]]),
      person_by_assay_variance_high = apply(boot_var, 2, function(z) sb_ci(z)[["high"]]),
      bootstrap_successful = colSums(is.finite(boot_pct)),
      bootstrap_attempted = B2_BOOT)
    if (boot_ci$bootstrap_successful[1] > 0L) {
      boot_ci$percent_explained_low[1] <- 0
      boot_ci$percent_explained_high[1] <- 0
    }
    if (any(boot_ci$bootstrap_successful < 0.90 * B2_BOOT))
      warning("B2 bootstrap success below 90%; inspect failed/undefined refits before reporting.")
    models <- models %>% left_join(boot_ci, by="model")
  }
  models <- models %>%
    mutate(ci_method = if (B2_BOOT > 0)
             "Participant-cluster percentile bootstrap (all models refitted per replicate)"
           else "No CI (point estimate only)",
           interpretation = paste("Exploratory descriptive variance reduction; anchor-matched modifiers;",
                                  "REML variance comparisons across fixed-effect specifications are not causal explained fractions"),
           CI_scope = "Conditional on the frozen calibration and complete-case modifier cohort")
} else {
  models <- tibble(model="not fitted",n_participants=n_distinct(d$RID),
                   reason="Fewer than 100 complete-case participants, inadequate repeat collections, or constant modifiers")
}
sb_save_csv(models,file.path(RES,"Table_S11_B2_variance_explanation.csv"),row.names=FALSE,na="")
cat("\nB2 coverage audit\n"); print(coverage)
cat("\nB2 complete-case variance models\n"); print(models)
.b2_products <- c("modifier_coverage", "modifier_value_audit", "matched_measurement_audit",
                  "unit_contract", "vitals_unit_audit", "variance_explanation",
                  if (model_identifiable && B2_BOOT > 0L) c("bootstrap_replicates", "bootstrap_diagnostics"))
for (.product in .b2_products)
  sb_stamp_result(file.path(RES, paste0("Table_S11_B2_", .product, ".csv")),
                  file.path(RAW, c("LABDATA.rda", "VITALS.rda")))
