if (!exists("sb_theme_classic",mode = "function"))
  stop("Run the updated 01_config.R first in this same R session.",call. = FALSE)
.sb_find <- function(file) {
  candidates <- c(if (exists("study_b_code_dir", inherits = TRUE))
    file.path(study_b_code_dir, file),
    file.path(Sys.getenv("STUDY_B_ROOT"), "CODE", file))
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in seq_len(7)) {
    candidates <- c(candidates, file.path(d, "CODE", file), file.path(d, file))
    if (identical(dirname(d), d)) break
    d <- dirname(d)
  }
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Cannot locate ", file, call. = FALSE)
  hit[1]
}
if (!exists("study_b_code_dir", inherits = TRUE))
  stop("Run 01_config.R manually first, in this same R session.", call. = FALSE)
if (!sb_have("sb_anchor", "sb_scored", "sb_rep_dat", "plasma_wide", "amy"))
  stop("Run 02_preprocess.R and 03_core_computation.R manually first.", call. = FALSE)
local({
  fail <- function(msg) stop("Figure S1 audit: ", msg, call. = FALSE)
  np <- function(d) dplyr::n_distinct(d$RID[!is.na(d$RID)])
  key <- c("RID", "plasma_date", "phase")
  ne <- function(d) nrow(dplyr::distinct(d, dplyr::across(dplyr::all_of(key))))
  check_equal <- function(a, b, columns, keys, label) {
    a <- a %>% select(all_of(columns)) %>% arrange(across(all_of(keys)))
    b <- b %>% select(all_of(columns)) %>% arrange(across(all_of(keys)))
    if (!isTRUE(all.equal(as.data.frame(a), as.data.frame(b),
                          check.attributes = FALSE, tolerance = 0)))
      fail(paste0(label, " differs from current analysis objects. Rebuild the ",
                  "pipeline in a clean session with current inputs/configuration; ",
                  "do not mix cached and newly reconstructed counts."))
  }
  raw <- sb_load_rda(file.path(study_b_data_dir,
                               "FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.rda"))
  required <- c("RID", "EXAMDATE", "ORIGPROT", "PLASMA_BIOMARKER", "TESTVALUE")
  if (length(setdiff(required, names(raw)))) fail("Missing required plasma columns.")
  selected <- raw %>% filter(PLASMA_BIOMARKER %in% unname(study_b_assay_map))
  mapped <- selected %>% transmute(
    RID = as.integer(as.character(RID)), plasma_date = sb_safe_date(EXAMDATE),
    phase = coalesce(na_if(toupper(trimws(as.character(ORIGPROT))), ""), "UNKNOWN"),
    assay = names(study_b_assay_map)[match(PLASMA_BIOMARKER, study_b_assay_map)],
    value = sb_safe_num(TESTVALUE, sentinel = FALSE),
    run_no = if ("RUN_NO" %in% names(raw)) sb_clean_id(RUN_NO) else NA_character_,
    sample_uid = if ("SAMPLEID" %in% names(raw)) sb_clean_id(SAMPLEID) else NA_character_,
    assay_cv = if ("CV" %in% names(raw)) sb_safe_num(CV, sentinel = FALSE) else NA_real_)
  id_ok <- mapped %>% filter(!is.na(RID))
  date_ok <- id_ok %>% filter(!is.na(plasma_date), !is.na(assay))
  dedup <- date_ok %>% distinct()
  conflicts <- dedup %>% group_by(across(all_of(c(key, "assay")))) %>%
    summarise(n_finite = n_distinct(value[is.finite(value)]), .groups = "drop") %>%
    filter(n_finite > 1)
  if (nrow(conflicts)) fail("Conflicting finite assay values; resolve in 02_preprocess.R.")
  collapsed <- dedup %>% group_by(across(all_of(c(key, "assay")))) %>%
    summarise(value = sb_first_finite(value), .groups = "drop")
  candidates <- collapsed %>% tidyr::pivot_wider(names_from = assay, values_from = value)
  for (a in setdiff(study_b_assays, names(candidates))) candidates[[a]] <- NA_real_
  finite <- candidates %>% filter(if_all(all_of(study_b_assays), ~is.finite(.x)))
  positive <- finite %>% filter(if_all(all_of(study_b_assays), ~.x > 0))
  check_equal(positive, plasma_wide, c(key, study_b_assays), key, "Plasma selection")
  check_equal(plasma_wide, sb_scored %>% distinct(RID, plasma_date, phase),
              key, key, "Scored collection set")
  stages <- list(mapped, id_ok, date_ok, dedup, collapsed, candidates, finite, positive)
  stage_names <- c("Selected four-assay source records", "Valid RID",
                   "Valid collection date and assay mapping", "Distinct transformed records",
                   "One value per RID/date/phase/assay", "Candidate collection occasions",
                   "All four assay values finite", "All four assay values strictly positive")
  ledger <- tibble(step = seq_along(stages), stage = stage_names,
                   unit = c(rep("assay records", 5), rep("collection occasions", 3)),
                   retained = vapply(stages, nrow, integer(1)),
                   participants = vapply(stages, np, integer(1))) %>%
    mutate(removed_from_previous = if_else(unit == lag(unit), lag(retained) - retained, NA_integer_),
           participants_lost_from_previous = lag(participants) - participants)
  reasons <- c("Start after selecting the four prespecified assays",
               "Missing or unparseable RID", "Missing/unparseable date or unmapped assay",
               "Identical transformed records collapsed (not participant exclusions)",
               "Repeated assay-cell records collapsed; conflicts stop processing",
               "Reshape: assay records to collection occasions (not an exclusion)",
               "At least one assay missing or non-finite",
               "At least one finite assay value zero or negative")
  ledger$reason <- reasons
  sb_write(ledger, "Figure_S1_plasma_selection_audit.csv")
  excluded_events <- candidates %>%
    mutate(exclusion_reason = case_when(
      !if_all(all_of(study_b_assays), ~is.finite(.x)) ~ "Missing/non-finite assay value",
      !if_all(all_of(study_b_assays), ~.x > 0) ~ "Non-positive finite assay value",
      TRUE ~ "Retained"))
  sb_write(excluded_events %>% select(all_of(key), exclusion_reason),
           "Figure_S1_collection_disposition.csv")
  pr <- sb_load_rda(file.path(study_b_data_dir, "UCBERKELEY_AMY_6MM.rda"))
  p0 <- pr %>% transmute(RID = as.integer(as.character(RID)), VISKEY = sb_viskey(VISCODE, VISCODE2),
                         date = sb_safe_date(SCANDATE),
                         amyloid_centiloid = sb_safe_num(CENTILOIDS, sentinel = FALSE),
                         amyloid_suvr = sb_safe_num(SUMMARY_SUVR, sentinel = FALSE),
                         amyloid_status = as.character(AMYLOID_STATUS),
                         qc_flag_raw = suppressWarnings(as.numeric(as.character(qc_flag))), qc_pass = sb_pet_qc_pass(qc_flag))
  p1 <- p0 %>% filter(!is.na(RID))
  p2 <- p1 %>% filter(!is.na(date))
  p3 <- p2 %>% filter(is.finite(amyloid_centiloid) | is.finite(amyloid_suvr))
  p4 <- p3 %>% distinct()
  p5 <- p4 %>% filter(qc_pass)
  p6 <- p5 %>% filter(is.finite(amyloid_centiloid))
  pc <- p6 %>% group_by(RID, date) %>% summarise(n = n_distinct(amyloid_centiloid), .groups="drop")
  if (any(pc$n > 1)) fail("Conflicting eligible PET Centiloids.")
  p7 <- p6 %>% arrange(RID, date, desc(qc_flag_raw)) %>% distinct(RID, date, .keep_all=TRUE)
  check_equal(p7, amy, c("RID", "date", "amyloid_centiloid"), c("RID", "date"), "PET selection")
  ps <- list(p0,p1,p2,p3,p4,p5,p6,p7)
  pet_ledger <- tibble(step=seq_along(ps), stage=c("Source PET records", "Valid RID",
                                                   "Valid scan date", "Finite Centiloid or SUVR", "Distinct transformed records",
                                                   "Allowed visual QC", "Finite Centiloid", "One record per RID/scan date"),
                       retained_records=vapply(ps,nrow,integer(1)), participants=vapply(ps,np,integer(1))) %>%
    mutate(removed_from_previous=lag(retained_records)-retained_records)
  sb_write(pet_ledger, "Figure_S1_PET_selection_audit.csv")
  reference_check <- sb_anchor %>%
    distinct(RID, sample_id, amyloid_centiloid, amyloid_positive) %>%
    mutate(expected_positive = if_else(is.finite(amyloid_centiloid),
                                       as.integer(amyloid_centiloid > study_b_centiloid_cutoff), NA_integer_))
  if (!identical(as.integer(reference_check$amyloid_positive),
                 as.integer(reference_check$expected_positive)))
    fail("Anchor PET classification differs from the configured Centiloid definition.")
  cutoff_label <- format(study_b_centiloid_cutoff, trim = TRUE, scientific = FALSE)
  anchor <- sb_anchor %>% distinct(RID, sample_id, amyloid_positive)
  if (anyDuplicated(anchor$RID)) fail("More than one primary anchor per participant.")
  rep_events <- sb_rep_dat %>% distinct(RID, sample_id)
  if (any((rep_events %>% count(RID))$n < 2)) fail("Longitudinal set contains singleton participants.")
  na <- nrow(anchor); nr <- np(rep_events); nv <- nrow(rep_events)
  nref <- sum(!is.na(anchor$amyloid_positive))
  npos <- sum(anchor$amyloid_positive == 1, na.rm=TRUE)
  nneg <- sum(anchor$amyloid_positive == 0, na.rm=TRUE)
  if (npos+nneg != nref) fail("PET classification must be binary or missing.")
  if (!study_b_anchor_requires_pet && na != np(positive)) fail("Anchor count mismatch.")
  numbers <- tibble(item=c("source_table_records", "selected_assay_records",
                           "candidate_occasions", "complete_occasions", "complete_participants",
                           "anchor_participants", "anchor_PET_available", "anchor_PET_unavailable",
                           "anchor_PET_positive", "anchor_PET_negative", "longitudinal_participants",
                           "longitudinal_total_occasions", "longitudinal_followup_occasions"),
                    n=c(nrow(raw),nrow(selected),nrow(candidates),nrow(positive),np(positive),
                        na,nref,na-nref,npos,nneg,nr,nv,nv-nr))
  sb_write(numbers,"Figure_S1_flowchart_numbers.csv")
  bx <- function(x,y,w,h,label) tibble(x=x,y=y,w=w,h=h,label=label)
  ed <- function(x,y,xend,yend) tibble(x=x,y=y,xend=xend,yend=yend)
  boxes <- bind_rows(
    bx(27,94,48,11,paste0("Supplied FNIH plasma data release\n",nrow(raw)," biomarker records")),
    bx(78,84,40,9,paste0("Other biomarker records not selected\n",nrow(raw)-nrow(selected)," records")),
    bx(27,74,48,13,paste0("Four prespecified assays\n",nrow(selected)," selected assay records\n",
                          nrow(candidates)," candidate collections; ",np(candidates)," participants")),
    bx(78,63,40,14,paste0("Collection exclusions\nMissing/non-finite: ",nrow(candidates)-nrow(finite),
                          "\nNon-positive finite: ",nrow(finite)-nrow(positive),
                          "\nParticipants lost: ",np(candidates)-np(positive))),
    bx(27,53,48,12,paste0("Complete four-assay cohort\nAll values finite and >0\n",
                          nrow(positive)," collections; ",np(positive)," participants")),
    bx(27,32,48,13,paste0("Primary classification analysis\nEarliest eligible plasma collection\n",
                          na," participants\n",if(study_b_anchor_requires_pet) "PET reference required" else "PET reference not required")),
    bx(78,32,40,13,paste0("Longitudinal persistence analysis\nAt least 2 eligible collections\n",
                          nr," participants; ",nv," collections\n",nv-nr," follow-up collections")),
    bx(78,19,40,8,paste0("Not included in PET performance\n",na-nref," anchors without eligible PET reference")),
    bx(27,8,48,12,paste0(
      "Amyloid PET-referenced analyses\n",
      nref," participants with eligible PET reference\n",
      npos," A-beta positive; ",nneg," A-beta negative\n",
      "Classification: Centiloid > ",cutoff_label
    )),
    bx(78,7,40,12,paste0(
      "PET reference matching\n",
      "Nearest eligible UCB amyloid PET scan\n",
      "Within +/-",study_b_amyloid_window," days\n",
      "Allowed QC and finite Centiloid\n",
      "Positive: Centiloid > ",cutoff_label
    )))
  edges <- bind_rows(ed(27,88.5,27,80.5),ed(27,84,58,84),ed(27,67.5,27,59),
                     ed(27,63,58,63),ed(27,47,27,38.5),ed(27,43,78,43),ed(78,43,78,38.5),
                     ed(27,25.5,27,14),ed(27,19,58,19))
  caption <- paste0(
    "Selection starts with the supplied release, not all ADNI recruits. ",
    "Collection = RID/date/phase. ",
    "Plasma records removed for invalid RID/date: ",
    nrow(mapped)-nrow(date_ok),
    "; duplicate-record reductions: ",
    nrow(date_ok)-nrow(collapsed),". ",
    "Analysis populations overlap. ",
    "Both binary and continuous PET analyses use the current QC-eligible, finite-Centiloid reference pool. ",
    "Binary positivity is defined as Centiloid > ",cutoff_label,
    "; AMYLOID_STATUS is not used for classification in this frozen analysis. ",
    "Anchors without an eligible reference are not necessarily participants without any PET scan. ",
    "Follow-up collection counts describe the repeated-collection cohort, not the final A3 comparison count. ",
    "PET reference pool before matching: ",
    nrow(p7),
    " records from ",
    np(p7),
    " participants. ",
    "QC-passing scans before finite-Centiloid filtering: ",
    nrow(p5),
    "; ",
    nrow(p5)-nrow(p6),
    " additionally removed for non-finite Centiloid. ",
    "Negative continuous Centiloids are retained."
  )
  writeLines(caption, file.path(study_b_dir, "Figure_S1_caption.txt"))
  figure_note <- paste0("Selection starts with the supplied FNIH release. Analysis populations overlap. ",
                        "PET reference requires allowed QC, finite Centiloid and matching within +/-",
                        study_b_amyloid_window," days; positivity: Centiloid > ",cutoff_label,
                        ". Missing eligible reference does not imply no PET scan. ",
                        "See Figure_S1_caption.txt and selection audit CSVs for detailed accounting.")
  p <- ggplot()+
    geom_segment(data=edges,aes(x=x,y=y,xend=xend,yend=yend),linewidth=.4,
                 arrow=grid::arrow(length=grid::unit(1.8,"mm")))+
    geom_rect(data=boxes,aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2),
              fill="white",colour="grey20",linewidth=.45)+
    geom_text(data=boxes,aes(x=x,y=y,label=label),size=3,lineheight=1.05)+
    coord_cartesian(xlim=c(0,100),ylim=c(0,102),expand=FALSE,clip="off")+
    labs(caption=paste(strwrap(figure_note,115),collapse="\n"))+
    theme_void()+theme(plot.caption=element_text(hjust=0,size=8),plot.margin=margin(12,12,12,12))
  sb_export_figure(p,file.path(study_b_dir,"Figure_S1_flowchart"),width = 190,height = 225)
  print(ledger,n=Inf);print(pet_ledger,n=Inf);print(numbers,n=Inf)
  message("Flow audit passed: reconstructed selections equal current analysis objects.")
})
