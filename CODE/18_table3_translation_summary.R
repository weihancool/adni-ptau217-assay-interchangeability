suppressPackageStartupMessages({ library(dplyr); library(purrr) })
root <- Sys.getenv("STUDY_B_ROOT")
res <- if (exists("study_b_dir")) study_b_dir else file.path(root, "RESULTS")
a1 <- sb_read_result(file.path(res, "Table_S7_A1_assay_switching_consequences.csv"), check.names = FALSE)
a2 <- sb_read_result(file.path(res, "Table_S8_A2_symptomatic_switching_consequences.csv"), check.names = FALSE)
a3 <- sb_read_result(file.path(res, "Table_S9_A3_longitudinal_pseudo_switch.csv"), check.names = FALSE)
tab <- bind_rows(
  a1 %>% filter(population == "Full cohort", scope == "Primary anchor", from_assay == "All") %>%
    transmute(section = "A1 Same-draw switching, full cohort", metric, n_participants,
              estimate_per_100 = per_100, ci_low = per_100_low, ci_high = per_100_high,
              estimand),
  a2 %>% filter(from_assay == "All") %>%
    transmute(section = "A2 Same-draw switching, MCI or dementia", metric, n_participants,
              estimate_per_100 = per_100, ci_low = per_100_low, ci_high = per_100_high,
              estimand = "Frozen full-cohort OOF decisions; equal-weight average over 12 directed switches"),
  a3 %>% filter(from_assay == "All", to_assay == "All",
                metric != "switch_alters_three_class_transition") %>%
    transmute(section = "A3 Longitudinal pseudo-switch", metric, n_participants,
              estimate_per_100 = per_100, ci_low = per_100_low, ci_high = per_100_high,
              estimand = case_when(
                metric == "same_assay_three_class_changed" ~
                  "Baseline assay A to follow-up assay A; three-class category",
                metric == "switched_three_class_changed" ~
                  "Baseline assay A to follow-up assay B; three-class category",
                metric == "category_change_excess_due_to_switch" ~
                  "Participant-equal switched minus same-assay three-class change rate",
                TRUE ~
                  "Baseline assay A to follow-up assay B on aligned same-draw first/last dates"
              ))
) %>%
  mutate(`Estimate (95% CI), per 100` = sprintf("%.1f (%.1f–%.1f)", estimate_per_100, ci_low, ci_high))
sb_save_csv(tab, file.path(res, "Table_3_translation_consequences.csv"), row.names = FALSE, na = "")
