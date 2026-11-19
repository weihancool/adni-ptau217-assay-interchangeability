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
if (!exists("sb_decision_long", inherits = TRUE))
  stop("Run 03_core_computation.R manually first.", call. = FALSE)
if (!exists("sb_vc_fig2_all", inherits=TRUE))
  stop("Run 09_figS5_variance_decomposition.R manually first.", call. = FALSE)
if (!exists("sb_variance_components_full", inherits=TRUE))
  stop("Run 10_tableS4_person_assay_audit.R manually first.", call. = FALSE)
sb_tableS4_components <- bind_rows(
  sb_vc_fig2_all,
  sb_variance_components_full %>%
    filter(!audit %in% c(
      "crossfit_main",
      "crossfit_plus_assay_by_fold"
    ))
)
stopifnot(
  !anyDuplicated(
    sb_tableS4_components[c("audit", "scale", "component")]
  )
)
stopifnot(!anyDuplicated(sb_tableS4_components[c("audit","scale","component")]))
sb_tableS4_components <- sb_tableS4_components %>%
  group_by(audit, scale) %>%
  group_modify(function(.x, .y) {
    get_variance <- function(pattern) {
      idx <- which(grepl(pattern, .x$component))
      if (length(idx) != 1L) {
        stop(
          "Variance component not uniquely identified: ",
          pattern, " | audit=", .y$audit,
          " | scale=", .y$scale,
          "\nAvailable components: ",
          paste(.x$component, collapse = ", "),
          call. = FALSE
        )
      }
      .x$variance[idx]
    }
    v_pa  <- get_variance("^person_by_assay$")
    v_occ <- get_variance("^occasion($|_)")
    v_res <- get_variance("^residual$")
    denom <- v_occ + v_res
    ratio <- if (
      all(is.finite(c(v_pa, v_occ, v_res))) &&
      all(c(v_pa, v_occ, v_res) >= 0) &&
      denom > 0
    ) {
      v_pa / denom
    } else {
      NA_real_
    }
    .x$person_assay_to_occasion_plus_residual_ratio <- ratio
    .x
  }) %>%
  ungroup()
sb_write(sb_tableS4_components,"Table_S4_variance_components_full.csv")
print(sb_tableS4_components,n=Inf)
