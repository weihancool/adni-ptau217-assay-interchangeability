suppressPackageStartupMessages({ library(dplyr); library(purrr); library(tidyr) })
root <- Sys.getenv("STUDY_B_ROOT")
res <- if (exists("study_b_dir")) study_b_dir else file.path(root, "RESULTS")
core <- sb_read_core(file.path(res, "02_core_results.rds"))
persons <- sb_read_result(file.path(res, "Table_1_participant_audit.csv"), check.names = FALSE)
symptomatic_ids <- persons %>%
  filter(diagnosis %in% c("MCI", "Dementia")) %>%
  transmute(RID = as.integer(RID)) %>% pull(RID)
d <- core$sb_anchor %>%
  mutate(RID = as.integer(as.character(RID)), assay = as.character(assay)) %>%
  filter(RID %in% symptomatic_ids, !is.na(amyloid_positive), !is.na(youden_positive))
exact_ci <- function(success, total) {
  if (!is.finite(total) || total <= 0) return(c(NA_real_, NA_real_))
  unname(stats::binom.test(success, total)$conf.int)
}
out <- d %>%
  group_by(assay) %>%
  summarise(
    total_n = n(),
    pet_positive_n = sum(amyloid_positive == 1L),
    pet_negative_n = sum(amyloid_positive == 0L),
    true_positive_n = sum(amyloid_positive == 1L & youden_positive == 1L),
    false_negative_n = sum(amyloid_positive == 1L & youden_positive == 0L),
    true_negative_n = sum(amyloid_positive == 0L & youden_positive == 0L),
    false_positive_n = sum(amyloid_positive == 0L & youden_positive == 1L),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    sensitivity = true_positive_n / pet_positive_n,
    sensitivity_low = exact_ci(true_positive_n, pet_positive_n)[1],
    sensitivity_high = exact_ci(true_positive_n, pet_positive_n)[2],
    specificity = true_negative_n / pet_negative_n,
    specificity_low = exact_ci(true_negative_n, pet_negative_n)[1],
    specificity_high = exact_ci(true_negative_n, pet_negative_n)[2],
    false_negative_rate = false_negative_n / pet_positive_n,
    false_negative_rate_low = exact_ci(false_negative_n, pet_positive_n)[1],
    false_negative_rate_high = exact_ci(false_negative_n, pet_positive_n)[2],
    false_positive_rate = false_positive_n / pet_negative_n,
    false_positive_rate_low = exact_ci(false_positive_n, pet_negative_n)[1],
    false_positive_rate_high = exact_ci(false_positive_n, pet_negative_n)[2],
    ci_method = "Clopper-Pearson exact binomial 95% CI",
    CI_scope = "Conditional on frozen OOF decisions; calibration uncertainty and shared training sets are not propagated"
  ) %>% ungroup()
frozen_path <- file.path(res, "Table_S8_A2_symptomatic_PET_performance.csv")
if (file.exists(frozen_path)) {
  frozen <- sb_read_result(frozen_path, check.names = FALSE)
  chk <- inner_join(
    out %>% select(assay, sensitivity, specificity, false_negative_rate, false_positive_rate),
    frozen %>% select(assay, sensitivity, specificity, false_negative_rate, false_positive_rate),
    by = "assay", suffix = c("_new", "_frozen")
  )
  cols <- c("sensitivity", "specificity", "false_negative_rate", "false_positive_rate")
  diffs <- unlist(lapply(cols, function(v) abs(chk[[paste0(v, "_new")]] - chk[[paste0(v, "_frozen")]])))
  stopifnot("A2 exact-CI comparison contains non-finite values" = all(is.finite(diffs)),
            "A2 exact-CI addition changed a frozen point estimate" =
              nrow(chk) == length(study_b_assays) && !anyDuplicated(chk$assay) &&
              length(diffs) > 0L && all(diffs < 1e-12))
}
sb_save_csv(out, file.path(res, "Table_S8_A2_PET_performance_exact_CI.csv"),
          row.names = FALSE, na = "")
message("Wrote A2 exact-binomial CI table.")
