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
if (!exists("sb_sensitivity", inherits = TRUE)) {
  stop("Run 21_figS4_sensitivity_window_cutoff.R manually first.", call. = FALSE)
}
if (!exists("sb_sensitivity", inherits = TRUE)) {
  stop(
    "26_tableS6: sensitivity results were not created by script 21.",
    call. = FALSE
  )
}
sb_write(
  sb_sensitivity,
  "Table_S6_sensitivity_window_cutoff_full.csv"
)
cat(
  "\nTable S6 -- window x cutoff sensitivity (full)\n"
)
print(as.data.frame(sb_sensitivity))
