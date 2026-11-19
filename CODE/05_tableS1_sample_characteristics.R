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
if (!sb_have("sb_anchor", "sb_rep_dat", "cohort"))
  stop("Run 02_preprocess.R and 03_core_computation.R manually first.", call. = FALSE)
if (!exists("table1_formatted", inherits=TRUE))
  stop("Run 04_table1_participant_characteristics.R manually first.", call. = FALSE)
sb_prevalence <- sb_anchor %>%
  pivot_longer(all_of(unname(sb_rules)), names_to = "rule", values_to = "positive") %>%
  mutate(rule = sub("_positive$", "", rule)) %>%
  group_by(rule, assay) %>%
  summarise(n = n(), positive_n = sum(positive),
            positive_percent = 100 * mean(positive), .groups = "drop")
sb_write(sb_prevalence, "Table_S1_assay_decision_prevalence.csv")
sb_category_prevalence <- sb_anchor %>%
  count(assay, three_class, name = "n") %>%
  group_by(assay) %>% mutate(percent = 100 * n / sum(n)) %>% ungroup()
sb_write(sb_category_prevalence, "Table_S1_three_class_prevalence.csv")
cat("\nTable S1 (lower half) -- binary decision prevalence\n"); print(sb_prevalence)
cat("\nTable S1 (lower half) -- three-class prevalence\n"); print(sb_category_prevalence)
if (exists("sb_tableS1_demographics", inherits = FALSE)) {
  cat("\nTable S1 -- demographics and time-matched diagnosis\n")
  print(sb_tableS1_demographics)
}
