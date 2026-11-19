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
if (!exists("sb_anchor", inherits = TRUE))
  stop("Run 03_core_computation.R manually first.", call. = FALSE)
.sb_fig1b_cache <- file.path(study_b_dir, "Figure_1B_pairwise_agreement_all_rules.csv")
sb_require_results("Figure_1B_pairwise_agreement_all_rules.csv",
                   "07_fig1_interchangeability.R",study_b_dir)
sb_pairwise <- sb_read_result(.sb_fig1b_cache,check.names = FALSE)
if (!all(c("rule","assay_1","assay_2","agreement","kappa","gwet_ac1") %in% names(sb_pairwise)) ||
    nrow(sb_pairwise) != length(sb_rules)*length(sb_pairs) ||
    anyDuplicated(sb_pairwise[c("rule","assay_1","assay_2")]))
  stop("Incompatible script 07 pairwise agreement output.",call. = FALSE)
sb_write(sb_pairwise, "Table_S2_pairwise_agreement_full.csv")
cat("\nTable S2 -- full pairwise agreement (18 rows)\n"); print(sb_pairwise)
