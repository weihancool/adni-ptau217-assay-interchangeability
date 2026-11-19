RESULTS_DIR <- Sys.getenv("STUDY_B_RESULTS_DIR", unset = file.path(Sys.getenv("STUDY_B_ROOT", unset = getwd()), "RESULTS"))
if (exists("study_b_dir", inherits = TRUE)) RESULTS_DIR <- get("study_b_dir", inherits = TRUE)
if (!dir.exists(RESULTS_DIR)) stop("RESULTS_DIR does not exist: ", RESULTS_DIR)
f <- file.path(RESULTS_DIR, "Table_S7_A1_residual_variance_benchmark.csv")
if (!file.exists(f)) stop("Missing: ", f, "\nRun updated 13 first.")
x <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
need <- c("source_margin_bin","comparisons","participants","observed_per100",
          "observed_low","observed_high","residual_null_per100",
          "excess_pp","excess_low","excess_high")
miss <- setdiff(need, names(x))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
fmt1 <- function(z) sprintf("%.1f", z)
fmt2small <- function(z) ifelse(abs(z) < 0.005, "<0.01", sprintf("%.2f", z))
fmt_ci1 <- function(e,l,h) sprintf("%.1f (%.1f–%.1f)", e,l,h)
panelD <- data.frame(
  `Absolute source margin, SD` = x$source_margin_bin,
  `People, n` = x$participants,
  `Contrasts, n` = format(x$comparisons, big.mark = ",", scientific = FALSE),
  `Observed binary reclassification per 100 (95% CI)` =
    fmt_ci1(x$observed_per100, x$observed_low, x$observed_high),
  `Residual-variance benchmark per 100` = fmt2small(x$residual_null_per100),
  `Excess over benchmark, percentage points (95% CI)` =
    fmt_ci1(x$excess_pp, x$excess_low, x$excess_high),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
csv_out <- file.path(RESULTS_DIR, "Table_S7_Panel_D_residual_variance_benchmark_display.csv")
tsv_out <- file.path(RESULTS_DIR, "Table_S7_Panel_D_residual_variance_benchmark_for_Word.tsv")
txt_out <- file.path(RESULTS_DIR, "Table_S7_Panel_D_residual_variance_benchmark_note.txt")
write.csv(panelD, csv_out, row.names = FALSE, na = "")
write.table(panelD, tsv_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
note <- paste0(
  "Panel D. Residual-variance benchmark for binary reclassification by source decision margin\n\n",
  "Note. Observed rates are the same all-occasion binary reclassification estimates reported in Panel C. ",
  "The residual-variance benchmark gives the expected threshold-crossing rate conditional on each observed source margin ",
  "if the destination-minus-source difference arose solely from two independent residual terms from the primary Youden-margin mixed model. ",
  "The residual term contains analytical error and other unmodeled assay-by-occasion variation and therefore should not be interpreted as pure analytical error. ",
  "Excess is observed minus benchmark; 95% confidence intervals were obtained from 2,000 participant-level cluster bootstrap replicates. ",
  "Participants may contribute comparisons to more than one margin stratum."
)
writeLines(note, txt_out, useBytes = TRUE)
print(panelD, row.names = FALSE)
cat("\nSaved:\n", csv_out, "\n", tsv_out, "\n", txt_out, "\n", sep = "")
