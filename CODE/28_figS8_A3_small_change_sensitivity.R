if (!exists("sb_theme_classic",mode = "function"))
  stop("Run the updated 01_config.R first in this same R session.",call. = FALSE)
suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(patchwork)})
res <- study_b_dir
sb_require_results(c("Table_S9_A3_small_change_threshold_sensitivity.csv",
  "Table_S9_A3_same_change_quartile_sensitivity.csv"),"15_A3_small_change_sensitivity.R",res)
thr <- sb_read_result(file.path(res,"Table_S9_A3_small_change_threshold_sensitivity.csv"),check.names = FALSE)
qrt <- sb_read_result(file.path(res,"Table_S9_A3_same_change_quartile_sensitivity.csv"),check.names = FALSE)
lab <- c(direction_reversal = "Direction reversal",platform_exceeds_same_change = "Platform > same-assay change")
cols <- c("Direction reversal" = sb_display_palette[["red"]],
          "Platform > same-assay change" = sb_display_palette[["blue"]])
for (z in list(thr,qrt)) {
  if (!all(c("metric","n_participants","per_100_estimate","per_100_low","per_100_high") %in% names(z)) ||
      any(!z$metric %in% names(lab))) stop("Incompatible A3 sensitivity inputs.",call. = FALSE)
}
keys <- sort(unique(thr$minimum_abs_same_change_SD))
n_by_cut <- vapply(keys,function(k) {
  n <- unique(thr$n_participants[thr$minimum_abs_same_change_SD == k])
  if (length(n) != 1L) stop("Metric denominators differ within a cutoff.",call. = FALSE)
  n
},numeric(1))
thr <- thr %>% mutate(label = lab[metric],
  stratum = factor(minimum_abs_same_change_SD,levels = rev(keys)))
qrt <- qrt %>% mutate(label = lab[metric],
  stratum = factor(same_change_quartile,levels = rev(sort(unique(same_change_quartile)))))
pd <- position_dodge(width = .32,orientation = "y")
make_plot <- function(z,letter,ylabel,ylabels) {
  ggplot(z,aes(per_100_estimate,stratum,colour = label,shape = label)) +
    geom_errorbar(aes(xmin = per_100_low,xmax = per_100_high),orientation = "y",
      width = .12,linewidth = .55,position = pd,na.rm = TRUE) +
    geom_point(size = 2,position = pd,na.rm = TRUE) +
    scale_colour_manual(values = cols) + scale_shape_manual(values = c(16,17)) +
    scale_x_continuous(limits = c(0,100),breaks = seq(0,100,25)) +
    scale_y_discrete(labels = ylabels,drop = FALSE) +
    labs(title = letter,x = "Events per 100 switches",y = ylabel,colour = NULL,shape = NULL) +
    sb_theme_classic() + theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(colour = sb_display_palette[["grid"]],linewidth = .3),
      legend.position = "bottom",legend.text = element_text(size = 7.5))
}
pA <- make_plot(thr,"A","Minimum |same-assay change| (SD)",
                setNames(sprintf("%.2f  (n = %d)",keys,n_by_cut),as.character(keys)))
pB <- make_plot(qrt,"B","Same-assay change quartile",
                function(x) paste0("Q",x,ifelse(x == "1"," (smallest)","")))
fig <- (pA | pB) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
sb_export_figure(fig,file.path(res,"Figure_S8_A3_small_change_sensitivity"),width = 190,height = 100)
sb_save_csv(bind_rows(thr %>% mutate(panel = "A"),qrt %>% mutate(panel = "B")),
            file.path(res,"Figure_S8_source_data.csv"),row.names = FALSE,na = "")
writeLines(paste0("Figure S8. Sensitivity of longitudinal counterfactual switching outcomes to same-assay ",
 "change magnitude. A: existing estimates/95% CIs under nested minimum-change restrictions; n is ",
 "participants retained, not switches. B: existing estimates/95% CIs by quartile, assigned once per ",
 "participant/source assay. Q1 contains the smallest absolute changes. The unrestricted row uses ",
 "the exact frozen A3 primary CI; restricted/quartile populations differ from the primary cohort. ",
 "Discrete strata are not interpolated; undefined strata are not zero rates."),file.path(res,"Figure_S8_caption.txt"))
if (interactive()) print(fig)
