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
if (!sb_have("sb_anchor"))
  stop("Run 03_core_computation.R manually first.", call. = FALSE)
sb_wide_class <- sb_anchor %>%
  select(RID, assay, three_class) %>%
  mutate(three_class = as.character(three_class)) %>%
  pivot_wider(names_from = assay, values_from = three_class)
stopifnot(!anyDuplicated(sb_wide_class$RID),
          all(study_b_assays %in% names(sb_wide_class)),
          !anyNA(sb_wide_class[, study_b_assays]))
.mat <- as.matrix(sb_wide_class[, study_b_assays])
sb_patterns <- sb_wide_class %>%
  mutate(pattern_type = case_when(
    apply(.mat, 1, function(z) length(unique(z)) == 1) ~ "complete agreement",
    apply(.mat, 1, function(z) "negative" %in% z && "positive" %in% z) ~
      "extreme negative-positive",
    TRUE ~ "gray-zone / adjacent only")) %>%
  unite("pattern", all_of(study_b_assays), sep = " | ", remove = FALSE) %>%
  count(pattern, pattern_type, sort = TRUE, name = "n") %>%
  mutate(percent = 100 * n / sum(n))
sb_write(sb_patterns, "Table_S3_gray_zone_patterns_full.csv")
sb_pattern_coverage <- tibble(
  observed_patterns = nrow(sb_patterns),
  theoretical_patterns = 3^length(study_b_assays),
  participants = nrow(sb_wide_class))
sb_write(sb_pattern_coverage, "Table_S3_pattern_coverage.csv")
sb_gray_pairwise <- map_dfr(sb_pairs, function(p) {
  x <- sb_wide_class[[p[1]]]; y <- sb_wide_class[[p[2]]]
  ind <- data.frame(agreement = as.numeric(x == y),
                    both_intermediate = as.numeric(x == "intermediate" &
                                                     y == "intermediate"),
                    either_intermediate = as.numeric(x == "intermediate" |
                                                       y == "intermediate"))
  n <- nrow(ind)
  set.seed(sb_seed(paste0("grayzone_", p[1], "_", p[2])))
  bo <- matrix(NA_real_, nrow = study_b_boot_n, ncol = 3,
               dimnames = list(NULL, names(ind)))
  for (b in seq_len(study_b_boot_n)) {
    i <- sample.int(n, n, replace = TRUE)
    bo[b, ] <- colMeans(ind[i, , drop = FALSE])
  }
  ci <- function(k) sb_ci(bo[, k])
  tibble(assay_1 = p[1], assay_2 = p[2], n = n,
         agreement = mean(ind$agreement),
         agreement_low = unname(ci("agreement")["low"]),
         agreement_high = unname(ci("agreement")["high"]),
         both_intermediate = mean(ind$both_intermediate),
         both_intermediate_low = unname(ci("both_intermediate")["low"]),
         both_intermediate_high = unname(ci("both_intermediate")["high"]),
         either_intermediate = mean(ind$either_intermediate),
         either_intermediate_low = unname(ci("either_intermediate")["low"]),
         either_intermediate_high = unname(ci("either_intermediate")["high"]),
         bootstrap_attempted = study_b_boot_n)
})
sb_write(sb_gray_pairwise, "Table_S3_gray_zone_pairwise.csv")
cat("\nTable S3 -- full three-class pattern table (", nrow(sb_patterns),
    " of ", 3^length(study_b_assays), " theoretical patterns observed)\n", sep = "")
print(sb_patterns, n = 60)
cat("\nTable S3 -- three-class pairwise agreement\n"); print(sb_gray_pairwise)
