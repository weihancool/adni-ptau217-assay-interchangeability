RUN_ORDER <- c(
  "01_config.R",
  "02_preprocess.R",
  "03_core_computation.R",
  "04_table1_participant_characteristics.R",
  "05_tableS1_sample_characteristics.R",
  "06_figS1_flowchart.R",
  "07_fig1_interchangeability.R",
  "08_table2_clinical_performance.R",
  "09_figS5_variance_decomposition.R",
  "10_tableS4_person_assay_audit.R",
  "11_tableS9_longitudinal_summary.R",
  "12_figS6_longitudinal_persistence.R",
  "13_clinical_switching_and_transport.R",
  "14_A2_pet_performance_exact_ci.R",
  "15_A3_small_change_sensitivity.R",
  "16_figure2_switching_consequences.R",
  "17_figure3_longitudinal_consequences.R",
  "18_table3_translation_summary.R",
  "19_figS2_gray_zone_patterns.R",
  "20_figS3_threshold_distributions.R",
  "21_figS4_sensitivity_window_cutoff.R",
  "22_tableS2_pairwise_agreement_full.R",
  "23_tableS4_variance_components_full.R",
  "24_tableS5_fixed_thresholds.R",
  "25_tableS3_gray_zone_patterns_full.R",
  "26_tableS6_sensitivity_full.R",
  "27_headline_synthesis.R",
  "28_figS8_A3_small_change_sensitivity.R",
  "29_B2_sample_key_audit.R",
  "30_B2_biological_modifiers.R",
  "31_final_freeze_checks.R"
)
print(data.frame(step = seq_along(RUN_ORDER),script = RUN_ORDER),row.names = FALSE)
if (identical(Sys.getenv("STUDY_B_RUN_ALL",unset = "0"),"1")) {
  root <- Sys.getenv("STUDY_B_ROOT",unset = "")
  if (!nzchar(root) || !file.exists(file.path(root,"CODE","01_config.R")))
    stop("Set STUDY_B_ROOT explicitly before an optional full run.",call. = FALSE)
  for (file in RUN_ORDER) {
    message("\n>>> RUNNING ",file)
    source(file.path(root,"CODE",file),local = .GlobalEnv,chdir = FALSE,encoding = "UTF-8")
  }
  message("Ordered analyses completed; inspect FINAL_FREEZE_CHECKS and FINAL_DISPLAY_EXPORT_CHECKS.")
} else message("Order reference only. Run 01-31 manually; no analysis was started.")
