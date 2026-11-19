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
if (!sb_have("sb_rep_dat", "sb_decision_long") ||
    !identical(get0("sb_core_contract_used", ifnotfound = NULL),
               sb_make_core_contract())) {
  stop("Run 03_core_computation.R manually first in this same R session.", call. = FALSE)
}
sb_spread_wide <- sb_rep_dat %>%
  distinct(RID, sample_id, plasma_date, assay_spread) %>%
  arrange(RID, plasma_date, sample_id) %>% group_by(RID) %>%
  summarise(first_spread = first(assay_spread),
            last_spread = last(assay_spread),
            interval_years = as.numeric(last(plasma_date) - first(plasma_date)) / 365.25,
            .groups = "drop")
sb_spread_rho <- suppressWarnings(
  cor(sb_spread_wide$first_spread, sb_spread_wide$last_spread, method = "spearman"))
set.seed(study_b_seed_spread)
sb_spread_boot <- replicate(study_b_boot_n, {
  i <- sample.int(nrow(sb_spread_wide), nrow(sb_spread_wide), replace = TRUE)
  suppressWarnings(cor(sb_spread_wide$first_spread[i], sb_spread_wide$last_spread[i],
                       method = "spearman"))
})
sb_spread_boot <- sb_spread_boot[is.finite(sb_spread_boot)]
sb_spread_ci <- sb_ci(sb_spread_boot)
sb_dev_wide <- sb_rep_dat %>%
  arrange(RID, assay, plasma_date, sample_id) %>% group_by(RID, assay) %>%
  filter(n() >= 2) %>%
  summarise(first_deviation = first(assay_deviation_adjusted),
            last_deviation = last(assay_deviation_adjusted),
            first_class = first(as.character(three_class)),
            last_class = last(as.character(three_class)), .groups = "drop")
sb_dev_rho_stat <- function(d) suppressWarnings(
  cor(d$first_deviation, d$last_deviation, method = "spearman", use = "complete.obs"))
sb_dev_sign_stat <- function(d) {
  ok <- !is.na(d$first_deviation) & !is.na(d$last_deviation)
  mean((d$first_deviation[ok] < 0) == (d$last_deviation[ok] < 0))
}
sb_dev_class_stat <- function(d) {
  ok <- !is.na(d$first_class) & !is.na(d$last_class)
  mean(d$first_class[ok] == d$last_class[ok])
}
sb_dev_rho_b <- sb_boot_summary(sb_dev_wide, sb_dev_rho_stat,
                                B = study_b_boot_n, seed = study_b_seed_dev)
sb_dev_sign_b <- sb_boot_summary(sb_dev_wide, sb_dev_sign_stat,
                                 B = study_b_boot_n, seed = study_b_seed_dev + 1L)
sb_dev_class_b <- sb_boot_summary(sb_dev_wide, sb_dev_class_stat,
                                  B = study_b_boot_n, seed = study_b_seed_dev + 2L)
sb_dev_by_assay <- map_dfr(study_b_assays, function(a) {
  d <- filter(sb_dev_wide, assay == a)
  r <- sb_boot_summary(d, sb_dev_rho_stat, B = study_b_boot_n, seed = study_b_seed_dev + 10L)
  s <- sb_boot_summary(d, sb_dev_sign_stat, B = study_b_boot_n, seed = study_b_seed_dev + 11L)
  cl <- sb_boot_summary(d, sb_dev_class_stat, B = study_b_boot_n,
                        seed = study_b_seed_dev + 12L)
  tibble(assay = a, n_participants = nrow(d),
         deviation_spearman = unname(r["estimate"]),
         deviation_spearman_low = unname(r["low"]),
         deviation_spearman_high = unname(r["high"]),
         sign_retention = unname(s["estimate"]),
         sign_retention_low = unname(s["low"]),
         sign_retention_high = unname(s["high"]),
         three_class_retention = unname(cl["estimate"]),
         three_class_retention_low = unname(cl["low"]),
         three_class_retention_high = unname(cl["high"]))
})
sb_write(sb_dev_by_assay, "Table_S9_supporting_deviation_by_assay.csv")
sb_rep_model_dat <- sb_rep_dat %>%
  mutate(phase = factor(phase), assay = factor(assay, levels = study_b_assays),
         RID = factor(RID))
sb_icc_stat <- function(d) {
  d <- d %>% mutate(phase = droplevels(factor(phase)),
                    assay = factor(assay, levels = study_b_assays))
  f <- assay_deviation ~ assay * years_from_baseline + (1|RID:assay)
  if (nlevels(d$phase) >= 2) f <- update(f, . ~ . + phase)
  m <- sb_checked_lmer(f,
            data = d, REML = TRUE,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 2e5)))
  v <- as.data.frame(VarCorr(m))
  vp <- v$vcov[v$grp == "RID:assay"][1]
  vr <- v$vcov[v$grp == "Residual"][1]
  sb_safe_ratio(vp, vp + vr)
}
sb_icc_point <- sb_icc_stat(sb_rep_model_dat)
sb_icc_ci <- c(low = NA_real_, high = NA_real_); sb_icc_nboot <- 0L
if (study_b_boot_n_model > 0) {
  message("Bootstrapping the deviation ICC with ", study_b_boot_n_model,
          " participant resamples.")
  sb_icc_boot_stat <- function(d) {
    id <- if (".boot_RID" %in% names(d)) ".boot_RID" else "RID"
    d2 <- d %>% mutate(RID = factor(.data[[id]]),
                       assay = factor(assay, levels = study_b_assays),
                       phase = droplevels(factor(phase)))
    sb_icc_stat(d2)
  }
  draws <- sb_cluster_boot(sb_rep_model_dat, sb_icc_boot_stat,
                           B = study_b_boot_n_model, seed = study_b_seed_dev + 20L)
  sb_icc_ci <- sb_ci(draws); sb_icc_nboot <- length(draws)
}
sb_persistence <- tibble(
  repeated_participants = length(unique(sb_rep_dat$RID)),
  median_interval_years = median(sb_spread_wide$interval_years),
  interval_q25 = quantile(sb_spread_wide$interval_years, .25),
  interval_q75 = quantile(sb_spread_wide$interval_years, .75),
  spread_spearman = sb_spread_rho,
  spread_low = unname(sb_spread_ci["low"]),
  spread_high = unname(sb_spread_ci["high"]),
  assay_deviation_spearman = unname(sb_dev_rho_b["estimate"]),
  assay_deviation_spearman_low = unname(sb_dev_rho_b["low"]),
  assay_deviation_spearman_high = unname(sb_dev_rho_b["high"]),
  assay_deviation_ICC = sb_icc_point,
  assay_deviation_ICC_low = unname(sb_icc_ci["low"]),
  assay_deviation_ICC_high = unname(sb_icc_ci["high"]),
  assay_deviation_ICC_bootstrap = sb_icc_nboot,
  assay_deviation_ICC_footnote = paste(
    "ICC from assay_deviation ~ assay*time + phase + (1|RID:assay);",
    "assay is fixed and the structurally inestimable occasion effect is omitted."),
  deviation_sign_retention = unname(sb_dev_sign_b["estimate"]),
  deviation_sign_retention_low = unname(sb_dev_sign_b["low"]),
  deviation_sign_retention_high = unname(sb_dev_sign_b["high"]),
  same_assay_three_class_retention = unname(sb_dev_class_b["estimate"]),
  same_assay_three_class_retention_low = unname(sb_dev_class_b["low"]),
  same_assay_three_class_retention_high = unname(sb_dev_class_b["high"]),
  bootstrap_successful = length(sb_spread_boot)
)
sb_write(sb_persistence, "Table_S9_supporting_longitudinal_persistence.csv")
cat("\nTable S9 supporting -- longitudinal persistence\n"); print(sb_persistence)
cat("\nTable S9 supporting -- per-assay deviation persistence\n"); print(sb_dev_by_assay)
local({
  ref <- sb_decision_long %>% filter(!is.na(amyloid_positive)) %>%
    arrange(RID, abs_amyloid_gap_days, plasma_date) %>% distinct(RID,assay,.keep_all=TRUE)
  constants <- map_dfr(study_b_assays,function(a) {
    r <- filter(ref,assay==a); mu <- mean(r$log_value); sig <- sd(r$log_value)
    th <- sb_pick_thresholds(r$amyloid_positive,(r$log_value-mu)/sig)
    tibble(assay=a,mu=mu,sig=sig,threshold=unname(th['youden']))
  })
  d <- sb_rep_dat %>% mutate(assay=as.character(assay)) %>%
    left_join(constants,by="assay",relationship="many-to-one") %>%
    mutate(fixed_margin=(log_value-mu)/sig-threshold) %>% group_by(RID,sample_id) %>%
    mutate(fixed_deviation=fixed_margin-mean(fixed_margin)) %>% ungroup()
  adjusted <- function(d, outcome, fold_adjust=FALSE) {
    terms <- "assay * years_from_baseline"
    if (n_distinct(d$phase)>1) terms <- c(terms,"factor(phase)")
    if (fold_adjust && n_distinct(d$fold)>1) terms <- c(terms,"assay:factor(fold)")
    f <- as.formula(paste(outcome,"~",paste(terms,collapse=" + ")))
    as.numeric(residuals(stats::lm(f,data=d,na.action=na.fail)))
  }
  d$fold_adjusted <- adjusted(d,"assay_deviation",TRUE)
  d$fixed_adjusted <- adjusted(d,"fixed_deviation",FALSE)
  modes <- c(crossfit_main="assay_deviation_adjusted",
             crossfit_assay_fold_adjusted="fold_adjusted",
             full_sample_fixed_calibration="fixed_adjusted")
  out <- imap_dfr(modes,function(col,mode) {
    w <- d %>% arrange(RID,assay,plasma_date,sample_id) %>% group_by(RID,assay) %>%
      filter(n()>=2) %>% summarise(first_deviation=first(.data[[col]]),
        last_deviation=last(.data[[col]]),.groups="drop")
    map_dfr(c("Pooled",study_b_assays),function(a) {
      z <- if(a=="Pooled") w else filter(w,assay==a)
      r <- sb_boot_summary(z,sb_dev_rho_stat,B=study_b_boot_n,
                           seed=sb_seed(paste0("persistence_audit_",a)))
      tibble(audit=mode,assay=a,participants=n_distinct(z$RID),
        rho=unname(r['estimate']),low=unname(r['low']),high=unname(r['high']),
        bootstrap_successful=unname(r['bootstrap_successful']),
        bootstrap_attempted=study_b_boot_n,
        CI_scope="Conditional on fitted calibration and deviation adjustment; fixed calibration is in-sample")
    })
  })
  sb_write(out,"Table_S9_supporting_persistence_calibration_audit.csv")
})
