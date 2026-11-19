suppressPackageStartupMessages(library(dplyr))
root <- Sys.getenv("STUDY_B_ROOT")
res <- if (exists("study_b_dir")) study_b_dir else file.path(root, "RESULTS")
core <- sb_read_core(file.path(res, "02_core_results.rds"))
events <- core$sb_events %>%
  transmute(RID = as.integer(as.character(RID)),
            sample_id = as.integer(as.character(sample_id)),
            plasma_date = as.Date(plasma_date))
audit <- tibble(
  check = c("sample_id globally unique in sb_events",
            "RID + sample_id unique in sb_events",
            "sample_id nonmissing", "RID nonmissing",
            "One collection per RID/date", "Collection dates nonmissing"),
  passed = c(!anyDuplicated(events$sample_id),
             !anyDuplicated(events[c("RID", "sample_id")]),
             !anyNA(events$sample_id), !anyNA(events$RID),
             !anyDuplicated(events[c("RID", "plasma_date")]), !anyNA(events$plasma_date)),
  n_violations = c(sum(duplicated(events$sample_id) | duplicated(events$sample_id, fromLast = TRUE)),
                   sum(duplicated(events[c("RID", "sample_id")]) |
                         duplicated(events[c("RID", "sample_id")], fromLast = TRUE)),
                   sum(is.na(events$sample_id)), sum(is.na(events$RID)),
                   sum(duplicated(events[c("RID", "plasma_date")])), sum(is.na(events$plasma_date)))
)
sb_save_csv(audit, file.path(res, "Table_S11_B2_sample_key_audit.csv"), row.names = FALSE)
stopifnot("B2 sample_id is not globally unique; existing joins are unsafe" = all(audit$passed))
message("B2 sample key audit passed.")
