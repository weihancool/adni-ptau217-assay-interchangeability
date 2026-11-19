if (!exists("sb_begin_reporting_step", mode = "function"))
  stop("Run revised 01_config.R first.", call. = FALSE)
sb_begin_reporting_step("32_reporting_input_audit.R")
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(lme4)})
res <- study_b_dir
core <- sb_read_core()
persons <- sb_read_result("Table_1_participant_audit.csv", check.names = FALSE)
s4 <- sb_read_result("Table_S4_variance_components_full.csv", check.names = FALSE)
fixed <- sb_read_result("Table_S4_fixed_threshold_audit.csv", check.names = FALSE)
d <- core$sb_decision_long
run_d <- d %>% group_by(occ_id) %>% filter(all(!is.na(run_no))) %>%
  ungroup() %>% mutate(assay_run = interaction(assay, run_no, drop = TRUE))
fixed_d <- d %>% mutate(assay_chr = as.character(assay)) %>%
  left_join(fixed, by = c("assay_chr" = "assay")) %>%
  mutate(fixed_global_z = (log_value - fixed_log_mean) / fixed_log_sd,
         fixed_margin_triage = fixed_global_z - fixed_triage_threshold_z,
         fixed_margin_confirmatory = fixed_global_z - fixed_confirmatory_threshold_z,
         fixed_margin_youden = fixed_global_z - fixed_youden_threshold_z)
first_pet <- d %>% arrange(RID, plasma_date, sample_id, assay) %>%
  group_by(RID) %>% summarise(first_pet = {
    z <- amyloid_positive[!is.na(amyloid_positive)]
    if (length(z)) first(z) else NA_integer_
  }, .groups = "drop") %>%
  mutate(baseline_amyloid_group = factor(case_when(
    first_pet == 1L ~ "positive", first_pet == 0L ~ "negative",
    TRUE ~ "unknown"), levels = c("negative","positive","unknown"))) %>%
  select(RID, baseline_amyloid_group)
first_pet_d <- d %>% select(-baseline_amyloid_group) %>%
  left_join(first_pet, by = "RID", relationship = "many-to-one")
keys <- s4 %>% distinct(audit, scale)
allowed <- c("crossfit_main","assay_random","crossfit_plus_assay_by_fold",
             "crossfit_run_subset_unadjusted","crossfit_plus_assay_by_run",
             "full_sample_fixed_threshold","first_available_PET_adjustment")
if (any(!keys$audit %in% allowed)) stop("Unknown S4 model specification.", call. = FALSE)
count_rows <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  audit <- keys$audit[i]; outcome <- keys$scale[i]
  dat <- switch(audit,
    crossfit_run_subset_unadjusted = run_d,
    crossfit_plus_assay_by_run = run_d,
    full_sample_fixed_threshold = fixed_d,
    first_available_PET_adjustment = first_pet_d,
    d)
  assay_fixed <- audit != "assay_random"
  terms <- c("(1|RID)","(1|occ_id)","(1|RID:assay)")
  if (!assay_fixed) terms <- c(terms,"(1|assay)")
  if (audit == "crossfit_plus_assay_by_fold") terms <- c(terms,"(1|assay_fold)")
  if (audit == "crossfit_plus_assay_by_run") terms <- c(terms,"(1|assay_run)")
  f <- as.formula(paste(outcome,"~",paste(c(sb_vc_fixed_terms(dat,assay_fixed),terms),collapse=" + ")))
  message("S4 sample-count validation ",i,"/",nrow(keys),": ",audit," / ",outcome)
  m <- sb_checked_lmer(f, data = dat, REML = TRUE,
    control = lme4::lmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5)))
  n <- sb_model_sample_counts(m)
  vc <- as.data.frame(lme4::VarCorr(m))
  group_map <- c(person="RID",occasion_within_person="occ_id",
    person_by_assay="RID:assay",assay="assay",assay_by_fold="assay_fold",
    assay_by_run="assay_run",residual="Residual")
  frozen <- s4 %>% filter(.data$audit == .env$audit, .data$scale == .env$outcome)
  frozen <- frozen[is.finite(frozen$variance),,drop=FALSE]
  fitted_v <- vapply(frozen$component,function(k) {
    z <- vc$vcov[vc$grp == group_map[[k]]]
    if (length(z) != 1L) stop("Variance component missing in validation refit: ",k)
    z
  },numeric(1))
  gap <- abs(fitted_v - frozen$variance)
  tolerance <- 1e-6 + 1e-5 * abs(frozen$variance)
  if (any(gap > tolerance))
    stop("Validation refit differs from frozen S4 estimates: ",audit," / ",outcome,
         ". No frozen estimate has been replaced.",call. = FALSE)
  for (k in intersect(c("participants","occasions","n_observations"),names(frozen))) {
    old <- frozen[[k]][is.finite(frozen[[k]])]
    if (length(old) && any(old != n[[k]])) stop("Frozen S4 count disagrees: ",audit," / ",outcome)
  }
  count_rows[[i]] <- data.frame(audit=audit,scale=outcome,
    participants=unname(n["participants"]),occasions=unname(n["occasions"]),
    n_observations=unname(n["n_observations"]),variance_match=TRUE,
    maximum_absolute_variance_difference=max(gap),
    count_source="model.frame/nobs from validation refit; frozen estimates and CIs retained")
}
sb_write(bind_rows(count_rows),"Table_S4_model_sample_counts_audit.csv")
pre_path <- file.path(res,"01_preprocessed.rds")
if (!identical(unname(tools::md5sum(pre_path)), .sb_frozen_preprocess_md5))
  stop("Preprocessing snapshot differs from the audited frozen release.",call. = FALSE)
.sb_note_input(pre_path)
pre <- readRDS(pre_path)
pet <- pre$amy %>% select(RID,date,amyloid_status)
if (anyDuplicated(pet[c("RID","date")])) stop("Frozen selected PET pool has duplicate RID/date.")
events <- core$sb_events %>%
  mutate(pet_selected_date=as.Date(plasma_date)+amyloid_gap_days) %>%
  left_join(pet,by=c("RID","pet_selected_date"="date"),suffix=c("","_pool"))
status <- if ("amyloid_status_pool" %in% names(events)) events$amyloid_status_pool else events$amyloid_status
status <- toupper(trimws(as.character(status)))
events$berkeley_positive <- case_when(
  status %in% c("1","POSITIVE","A+","TRUE") ~ 1L,
  status %in% c("0","NEGATIVE","A-","FALSE") ~ 0L, TRUE ~ NA_integer_)
concordance <- bind_rows(
  events %>% mutate(scope="All matched collection occasions"),
  events %>% filter(primary_anchor) %>% mutate(scope="Primary anchors")) %>%
  filter(!is.na(amyloid_positive)) %>%
  count(scope,CL_gt20=amyloid_positive,AMYLOID_STATUS=berkeley_positive,name="n")
sb_write(concordance,"FINAL_PET_CL20_STATUS_CROSSTAB.csv")
sb_write(events %>% filter(!is.na(amyloid_positive),!is.na(berkeley_positive),
  amyloid_positive != berkeley_positive) %>%
  select(RID,sample_id,primary_anchor,plasma_date,pet_selected_date,
         amyloid_centiloid,amyloid_positive,berkeley_positive),
  "FINAL_PET_CL20_STATUS_DISCORDANCE.csv")
window_counts <- bind_rows(lapply(c(90,180,365),function(w) data.frame(
  window_days=w,anchor_participants=nrow(persons),
  diagnosis_available=sum(is.finite(persons$diagnosis_gap_days)&persons$diagnosis_gap_days<=w),
  MCI_or_dementia=sum(persons$diagnosis %in% c("MCI","Dementia") &
                      is.finite(persons$diagnosis_gap_days)&persons$diagnosis_gap_days<=w),
  MMSE_available=sum(is.finite(persons$mmse_gap_days)&persons$mmse_gap_days<=w),
  CDRSB_available=sum(is.finite(persons$cdrsb_gap_days)&persons$cdrsb_gap_days<=w))))
sb_write(window_counts,"FINAL_CLINICAL_MATCH_WINDOW_AUDIT.csv")
sb_write(persons %>% filter(diagnosis_gap_days>180 | mmse_gap_days>180 | cdrsb_gap_days>180) %>%
  select(RID,diagnosis,diagnosis_gap_days,mmse_gap_days,cdrsb_gap_days),
  "FINAL_CLINICAL_MATCH_OVER180_AUDIT.csv")
read_raw <- function(name) {
  path <- file.path(study_b_data_dir,paste0(name,".rda"))
  if (!file.exists(path)) stop("Raw audit input absent: ",path,call. = FALSE)
  .sb_note_input(path); sb_load_rda(path,name)
}
fnih <- read_raw("FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES")
raw_values <- fnih %>% transmute(RID=as.integer(as.character(RID)),
  plasma_date=sb_safe_date(EXAMDATE),
  enrollment_protocol=as.character(ORIGPROT),
  assay=names(study_b_assay_map)[match(PLASMA_BIOMARKER,study_b_assay_map)],
  reported_value=as.character(TESTVALUE),
  numeric_value=sb_safe_num(TESTVALUE,sentinel=FALSE)) %>%
  filter(!is.na(assay),!is.na(RID),!is.na(plasma_date)) %>%
  mutate(boundary_string=grepl("[<>]|LLOQ|ULOQ|BLQ|BELOW|ABOVE",toupper(reported_value)),
         finite_positive=is.finite(numeric_value)&numeric_value>0)
sb_write(raw_values %>% group_by(assay,reported_value,boundary_string,finite_positive) %>%
  summarise(source_rows=n(),.groups="drop"),"FINAL_FNIH_VALUE_ENCODING_AUDIT.csv")
wide <- raw_values %>% group_by(RID,plasma_date,enrollment_protocol,assay) %>%
  summarise(value=sb_first_finite(numeric_value),.groups="drop") %>%
  pivot_wider(names_from=assay,values_from=value)
excluded <- wide %>% filter(!if_all(all_of(study_b_assays),~is.finite(.x)&.x>0))
sb_write(excluded,"FINAL_FNIH_EXCLUDED_COLLECTION_AUDIT.csv")
flag_names <- names(fnih)[grepl("LLOQ|ULOQ|LOQ|QC|COMMENT|NOTE|FLAG",names(fnih),ignore.case=TRUE)]
sb_write(data.frame(field=flag_names),"FINAL_FNIH_RANGE_FLAG_FIELDS.csv")
if (length(flag_names)) sb_write(fnih %>%
  filter(PLASMA_BIOMARKER %in% unname(study_b_assay_map)) %>%
  select(any_of(c("RID","EXAMDATE","PLASMA_BIOMARKER","TESTVALUE",flag_names))),
  "FINAL_FNIH_RANGE_FLAG_RECORDS.csv")
amy_raw <- read_raw("UCBERKELEY_AMY_6MM")
if ("TRACER" %in% names(amy_raw)) {
  tracers <- amy_raw %>% transmute(RID=as.integer(as.character(RID)),
    pet_selected_date=sb_safe_date(SCANDATE),TRACER=as.character(TRACER)) %>% distinct()
  tracer_audit <- events %>% filter(!is.na(amyloid_positive)) %>%
    select(RID,sample_id,primary_anchor,pet_selected_date) %>%
    left_join(tracers,by=c("RID","pet_selected_date"),relationship="many-to-many")
  sb_write(tracer_audit,"FINAL_SELECTED_PET_TRACER_AUDIT.csv")
} else sb_write(data.frame(status="TRACER field absent; no tracer was inferred"),
               "FINAL_SELECTED_PET_TRACER_AUDIT.csv")
vital_path <- file.path(study_b_data_dir,"VITALS.rda")
if (file.exists(vital_path)) {
  vitals <- read_raw("VITALS")
  unit_rows <- bind_rows(lapply(intersect(c("VSWTUNIT","VSHTUNIT"),names(vitals)),function(k)
    data.frame(field=k,raw_unit=as.character(vitals[[k]])))) %>%
    count(field,raw_unit,name="source_rows") %>%
    mutate(requires_release_dictionary=grepl("^[0-9]+$",raw_unit))
  sb_write(unit_rows,"FINAL_VITALS_UNIT_ENCODING_AUDIT.csv")
} else sb_write(data.frame(status="VITALS input absent; units not inferred"),
               "FINAL_VITALS_UNIT_ENCODING_AUDIT.csv")
sb_write(data.frame(policy="Exact audited frozen release; original sidecars retained",
  frozen_core_md5=.sb_frozen_core_md5,
  current_config_md5=unname(tools::md5sum(file.path(study_b_code_dir,"01_config.R"))),
  original_config_md5=.sb_frozen_code_md5[["01_config.R"]],
  original_analysis_identity_preserved=TRUE),"FINAL_PROVENANCE_COMPATIBILITY.csv")
sb_end_reporting_step()
message("Reporting input audit completed. Frozen statistical estimates/CIs were retained.")
