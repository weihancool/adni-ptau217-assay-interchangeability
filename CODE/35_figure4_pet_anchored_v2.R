suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})
ASSAYS <- c("AlzPath", "C2N_ratio", "Fujirebio", "Janssen")
DISPLAY <- c(AlzPath="ALZpath", C2N_ratio="C2N", Fujirebio="Fujirebio", Janssen="Janssen")
SOURCE_BD <- "C2N_ratio"
PET_WINDOW <- 90L
PET_CUTOFF <- 20
BOOT_N <- 1000L
BOOT_SEED <- 20261016L
FIG_WIDTH_MM <- 190
FIG_HEIGHT_MM <- 142
SAVE_DPI <- 600
ink   <- "#253444"
blue  <- "#254F70"
red   <- "#B44145"
pale  <- "#E7EBEE"
grey  <- "#B2BCC5"
white <- "#FFFFFF"
fig_theme <- theme_classic(base_size=9) +
  theme(
    text=element_text(colour=ink),
    axis.title=element_text(colour=ink, size=8),
    axis.text=element_text(colour=ink, size=7.5),
    axis.line=element_line(colour=ink, linewidth=.35),
    axis.ticks=element_line(colour=ink, linewidth=.35),
    plot.title=element_text(face="bold", size=12, hjust=0, margin=margin(b=2)),
    plot.title.position="plot",
    plot.subtitle=element_blank(),
    strip.background=element_blank(),
    strip.text=element_text(colour=ink, face="bold", size=7.5),
    legend.key=element_blank(),
    legend.text=element_text(size=7),
    plot.margin=margin(3,5,3,4)
  )
find_results_dir <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  candidates <- c(
    Sys.getenv("STUDY_B_RESULTS_DIR", unset=""),
    if (exists("study_b_dir", inherits=TRUE)) get("study_b_dir", inherits=TRUE) else character(),
    if (nzchar(Sys.getenv("STUDY_B_ROOT", unset="")))
      file.path(Sys.getenv("STUDY_B_ROOT"), "RESULTS") else character(),
    if (length(args)) args[[1]] else character(),
    file.path(getwd(), "RESULTS"), file.path(dirname(getwd()), "RESULTS")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  found <- candidates[file.exists(file.path(candidates, "02_core_results.rds"))]
  if (!length(found)) stop(
    "Cannot locate RESULTS/02_core_results.rds. Set STUDY_B_RESULTS_DIR.\nChecked: ",
    paste(candidates, collapse=" | "), call.=FALSE)
  normalizePath(found[[1]], winslash="/", mustWork=TRUE)
}
res <- find_results_dir()
out <- file.path(res, "PET_FIG4")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write_csv <- function(x, filename) {
  utils::write.csv(x, file.path(out, filename), row.names=FALSE, na="")
}
core <- readRDS(file.path(res, "02_core_results.rds"))
if (!is.list(core) || is.null(core$sb_scored))
  stop("02_core_results.rds lacks core$sb_scored.", call.=FALSE)
d <- as_tibble(core$sb_scored)
required <- c("RID","sample_id","plasma_date","assay","amyloid_positive",
              "amyloid_centiloid","abs_amyloid_gap_days","youden_positive")
miss <- setdiff(required, names(d))
if (length(miss)) stop("Missing core columns: ", paste(miss, collapse=", "), call.=FALSE)
d <- d %>% mutate(
  assay=as.character(assay), plasma_date=as.Date(plasma_date),
  amyloid_positive=as.integer(amyloid_positive),
  amyloid_centiloid=as.numeric(amyloid_centiloid),
  abs_amyloid_gap_days=as.numeric(abs_amyloid_gap_days),
  youden_positive=as.integer(youden_positive))
if (!setequal(unique(d$assay), ASSAYS)) stop("Unexpected assay names.", call.=FALSE)
if (anyNA(d$RID) || anyNA(d$sample_id) || anyNA(d$plasma_date))
  stop("Missing IDs or dates.", call.=FALSE)
if (anyDuplicated(d[c("RID","sample_id","assay")]))
  stop("Duplicate RID/sample_id/assay rows.",call.=FALSE)
if (any(table(d$sample_id) != 4L))
  stop("Some blood draws lack all four assays.", call.=FALSE)
if (anyNA(d$youden_positive) || any(!(d$youden_positive %in% 0:1)))
  stop("Invalid frozen binary classifications.", call.=FALSE)
pet_rows <- !is.na(d$amyloid_positive)
if (any(!d$amyloid_positive[pet_rows] %in% 0:1) ||
    any(!is.finite(d$amyloid_centiloid[pet_rows])) ||
    any(!is.finite(d$abs_amyloid_gap_days[pet_rows])))
  stop("Invalid PET reference/gap/Centiloid.", call.=FALSE)
if (any(d$amyloid_positive[pet_rows] !=
        as.integer(d$amyloid_centiloid[pet_rows] > PET_CUTOFF)))
  stop("Core PET status differs from Centiloid > 20.", call.=FALSE)
fl <- d %>% arrange(RID, assay, plasma_date, sample_id) %>%
  group_by(RID, assay) %>% summarise(
    n_draws=n(),
    first_sample=first(sample_id), last_sample=last(sample_id),
    first_date=first(plasma_date), last_date=last(plasma_date),
    first_pet=first(amyloid_positive), last_pet=last(amyloid_positive),
    first_cl=first(amyloid_centiloid), last_cl=last(amyloid_centiloid),
    first_gap=first(abs_amyloid_gap_days), last_gap=last(abs_amyloid_gap_days),
    first_binary=first(youden_positive), last_binary=last(youden_positive),
    .groups="drop")
audit <- fl %>% group_by(RID) %>% summarise(
  n_assays=n_distinct(assay), n_min=min(n_draws), n_max=max(n_draws),
  n_first_sample=n_distinct(first_sample, na.rm=FALSE),
  n_last_sample=n_distinct(last_sample, na.rm=FALSE),
  n_first_date=n_distinct(first_date, na.rm=FALSE),
  n_last_date=n_distinct(last_date, na.rm=FALSE),
  n_first_pet=n_distinct(first_pet, na.rm=FALSE),
  n_last_pet=n_distinct(last_pet, na.rm=FALSE),
  n_first_cl=n_distinct(first_cl, na.rm=FALSE),
  n_last_cl=n_distinct(last_cl, na.rm=FALSE),
  n_first_gap=n_distinct(first_gap, na.rm=FALSE),
  n_last_gap=n_distinct(last_gap, na.rm=FALSE), .groups="drop")
bad <- audit %>% filter(n_assays!=4 | n_min<2 | n_min!=n_max |
  n_first_sample!=1 | n_last_sample!=1 | n_first_date!=1 | n_last_date!=1 |
  n_first_pet!=1 | n_last_pet!=1 | n_first_cl!=1 | n_last_cl!=1 |
  n_first_gap!=1 | n_last_gap!=1)
if (nrow(bad)) {
  write_csv(bad, "Figure_4_AUDIT_invalid_panels.csv")
  stop("Assay panels not aligned; see Figure_4_AUDIT_invalid_panels.csv.", call.=FALSE)
}
persons <- fl %>% filter(assay == ASSAYS[1]) %>%
  filter(!is.na(first_pet), !is.na(last_pet),
         is.finite(first_gap), is.finite(last_gap),
         first_gap <= PET_WINDOW, last_gap <= PET_WINDOW,
         is.finite(first_cl), is.finite(last_cl)) %>%
  mutate(
    PET_transition=case_when(
      first_pet==0 & last_pet==0 ~ "PET 0->0",
      first_pet==0 & last_pet==1 ~ "PET 0->1",
      first_pet==1 & last_pet==1 ~ "PET 1->1",
      first_pet==1 & last_pet==0 ~ "PET 1->0"),
    years_between=as.numeric(last_date-first_date)/365.25)
if (!nrow(persons)) stop("No participants have PET at both ends.", call.=FALSE)
if (anyNA(persons$PET_transition)) stop("Unexpected PET state.", call.=FALSE)
coverage <- persons %>% count(PET_transition, name="n")
coverage <- tibble(PET_transition=c("PET 0->0","PET 0->1","PET 1->1","PET 1->0")) %>%
  left_join(coverage,by="PET_transition") %>% mutate(n=replace_na(n,0L))
write_csv(coverage,"Figure_4A_PET_transitions.csv")
pairs <- map_dfr(ASSAYS, function(src) map_dfr(setdiff(ASSAYS,src), function(dst) {
  a <- fl %>% filter(assay==src)
  b <- fl %>% filter(assay==dst)
  ab <- inner_join(a,b,by="RID",suffix=c("_A","_B"))
  if (nrow(ab)!=nrow(a) ||
      any(ab$first_sample_A != ab$first_sample_B) ||
      any(ab$last_sample_A != ab$last_sample_B) ||
      any(ab$first_date_A != ab$first_date_B) ||
      any(ab$last_date_A != ab$last_date_B))
    stop("First/last samples misaligned for ",src,"->",dst,call.=FALSE)
  ab %>% transmute(
    RID, source=src, destination=dst,
    first_pet=first_pet_A,last_pet=last_pet_A,
    first_cl=first_cl_A,last_cl=last_cl_A,
    first_gap=first_gap_A,last_gap=last_gap_A,
    A_first=first_binary_A,A_last=last_binary_A,B_last=last_binary_B)
})) %>%
  semi_join(persons %>% select(RID),by="RID") %>%
  mutate(
    source_first_correct=A_first==first_pet,
    new_discordance=(A_last==last_pet & B_last!=last_pet),
    resolved_discordance=(A_last!=last_pet & B_last==last_pet),
    net_plasma_positive=as.integer(B_last==1)-as.integer(A_last==1),
    classification_status=case_when(
      new_discordance ~ "New PET discordance",
      resolved_discordance ~ "Resolved PET discordance",
      TRUE ~ "Unchanged PET concordance"))
if (nrow(pairs)!=12L*nrow(persons) || anyDuplicated(pairs[c("RID","source","destination")]))
  stop("Expected 12 unique directed comparisons per eligible participant.",call.=FALSE)
stable <- pairs %>% filter(first_pet==0L,last_pet==0L,
                           source_first_correct, A_first==0L)
expected_stable <- coverage$n[coverage$PET_transition=="PET 0->0"]
if (!nrow(stable) || n_distinct(stable$RID) > expected_stable)
  stop("Invalid PET-stable negative subgroup.",call.=FALSE)
pair_flow <- stable %>% group_by(source,destination) %>% summarise(
  n=n(),
  agree_negative=sum(A_last==0L & B_last==0L),
  new_positive=sum(A_last==0L & B_last==1L),
  resolved_positive=sum(A_last==1L & B_last==0L),
  agree_positive=sum(A_last==1L & B_last==1L),
  positive_keep=sum(A_last==1L), positive_switch=sum(B_last==1L),
  .groups="drop") %>% mutate(
    net_change_pp=100*(new_positive-resolved_positive)/n,
    keep_positive_percent=100*positive_keep/n,
    switch_positive_percent=100*positive_switch/n)
if (nrow(pair_flow)!=12L || any(pair_flow$n!=pair_flow$agree_negative+
    pair_flow$new_positive+pair_flow$resolved_positive+pair_flow$agree_positive))
  stop("Pairwise 2x2 counts failed verification.",call.=FALSE)
write_csv(pair_flow,"Figure_4BC_stable_negative_pair_flows.csv")
source_person <- stable %>% group_by(source,RID) %>% summarise(
  person_net=mean(net_plasma_positive), .groups="drop")
ci_cluster <- function(x, seed) {
  x <- as.numeric(x)
  n <- length(x)
  est <- 100*mean(x)
  if (!n || BOOT_N < 1L || n < 10L)
    return(c(estimate=est,low=NA_real_,high=NA_real_))
  set.seed(seed)
  counts <- stats::rmultinom(BOOT_N, size=n, prob=rep(1/n,n))
  draws <- 100*as.numeric(crossprod(x,counts))/n
  c(estimate=est,
    low=unname(quantile(draws,.025,names=FALSE)),
    high=unname(quantile(draws,.975,names=FALSE)))
}
source_ci <- map_dfr(seq_along(ASSAYS),function(i) {
  a <- ASSAYS[i]
  x <- source_person %>% filter(source==a)
  z <- ci_cluster(x$person_net, BOOT_SEED+i-1L)
  tibble(source=a, n_participants=nrow(x),
         mean_net_pp=unname(z["estimate"]), ci_low=unname(z["low"]),
         ci_high=unname(z["high"]))
})
write_csv(source_ci,"Figure_4BC_source_averages_bootstrap.csv")
for (a in ASSAYS) {
  v <- pair_flow %>% filter(source==a)
  z <- source_ci %>% filter(source==a)
  if (abs(mean(v$net_change_pp)-z$mean_net_pp)>1e-8)
    stop("Source average does not match equal-weight pair average for ",a)
  if (any(v$n!=z$n_participants))
    stop("One or more destinations have a different denominator for ",a)
}
stable_far <- stable %>% filter(abs(first_cl-PET_CUTOFF)>5,
                               abs(last_cl-PET_CUTOFF)>5)
far_stats <- stable_far %>% group_by(source,RID) %>%
  summarise(person_net=mean(net_plasma_positive),.groups="drop") %>%
  group_by(source) %>% summarise(n_participants=n(),
                                mean_net_pp=100*mean(person_net),.groups="drop")
write_csv(far_stats,"Figure_4_sensitivity_exclude_CL15_to25.csv")
transition_n <- setNames(coverage$n,coverage$PET_transition)
median_years <- median(persons$years_between)
boxes <- tibble(
  xmin=c(.30,5.25,.30,5.25), xmax=c(4.72,9.67,4.72,9.67),
  ymin=c(4.60,4.60,3.27,3.27), ymax=c(5.45,5.45,4.12,4.12),
  label=c("First PET","Last PET","First blood: assay A","Last blood: A or B"))
state_boxes <- tibble(
  xmin=c(.30,5.25,.30,5.25), xmax=c(4.72,9.67,4.72,9.67),
  ymin=c(1.20,1.20,.15,.15), ymax=c(2.09,2.09,1.04,1.04),
  label=c(sprintf("PET - / -     n = %d",transition_n[["PET 0->0"]]),
          sprintf("PET - / +     n = %d",transition_n[["PET 0->1"]]),
          sprintf("PET + / +     n = %d",transition_n[["PET 1->1"]]),
          sprintf("PET + / -     n = %d",transition_n[["PET 1->0"]])))
pA <- ggplot() +
  geom_rect(data=boxes,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill="#F1F4F6",colour="#A7B2BC",linewidth=.35) +
  geom_text(data=boxes,aes(x=(xmin+xmax)/2,y=(ymin+ymax)/2,label=label),
            size=2.75,colour=ink) +
  annotate("segment",x=4.79,xend=5.17,y=5.02,yend=5.02,
           arrow=grid::arrow(length=grid::unit(1.5,"mm")),linewidth=.4,colour=ink) +
  annotate("segment",x=4.79,xend=5.17,y=3.69,yend=3.69,
           arrow=grid::arrow(length=grid::unit(1.5,"mm")),linewidth=.4,colour=ink) +
  annotate("text",x=5,y=5.83,label=sprintf("Median %.1f years",median_years),
           size=2.45,colour=ink) +
  annotate("text",x=.33,y=2.59,label=sprintf("Paired PET: n = %d",nrow(persons)),
           size=2.6,hjust=0,colour=ink) +
  annotate("segment",x=.30,xend=9.67,y=2.39,yend=2.39,
           colour=pale,linewidth=.4) +
  geom_rect(data=state_boxes,
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill="#F7F9FA",colour="#D4DCE2",linewidth=.3) +
  geom_text(data=state_boxes,
            aes(x=(xmin+xmax)/2,y=(ymin+ymax)/2,label=label),
            size=2.6,colour=ink) +
  coord_cartesian(xlim=c(0,10),ylim=c(0,6.15),expand=FALSE,clip="off") +
  labs(title="A") + fig_theme +
  theme(axis.title=element_blank(),axis.text=element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),
        plot.margin=margin(3,8,2,5))
focus <- stable %>% filter(source==SOURCE_BD)
focus_n <- n_distinct(focus$RID)
if (nrow(focus)!=3L*focus_n || focus_n < 1L)
  stop("B: C2N participants lack three destination classifications.",call.=FALSE)
focus_avg <- source_ci %>% filter(source==SOURCE_BD)
dest <- setdiff(ASSAYS,SOURCE_BD)
focus_flows <- pair_flow %>% filter(source==SOURCE_BD)
facet_map <- setNames(unname(DISPLAY[focus_flows$destination]),
                      focus_flows$destination)
facet_levels <- unname(facet_map[dest])
matrix_points <- focus %>% count(destination,A_last,B_last,name="n")
mx <- tidyr::expand_grid(destination=dest,A_last=0:1,B_last=0:1) %>%
  left_join(matrix_points,by=c("destination","A_last","B_last")) %>%
  mutate(n=replace_na(n,0L),
         transition=case_when(
           A_last==0L & B_last==1L ~ "New PET discordance",
           A_last==1L & B_last==0L ~ "Resolved PET discordance",
           TRUE ~ "Unchanged classification"),
         row=factor(if_else(A_last==0L,"-","+"),levels=c("+","-")),
         col=factor(if_else(B_last==0L,"-","+"),levels=c("-","+")),
         facet=factor(facet_map[destination],levels=facet_levels))
if (any(mx$n < 0) || anyNA(mx$n) || nrow(mx)!=12L ||
    any(tapply(mx$n,mx$destination,sum)!=focus_n))
  stop("B: paired-matrix cell counts are inconsistent.",call.=FALSE)
write_csv(mx %>% mutate(row=as.character(row),col=as.character(col),
                         facet=as.character(facet)),"Figure_4B_C2N_paired_matrices.csv")
fill_class <- c("New PET discordance"=red,
                "Resolved PET discordance"=blue,
                "Unchanged classification"="#E2E7EB")
pB <- ggplot(mx,aes(col,row,fill=transition)) +
  geom_tile(colour=white,linewidth=1.2,width=.98,height=.98) +
  geom_text(aes(label=n, colour=transition),size=3.15,
            fontface="plain",show.legend=FALSE) +
  facet_wrap(~facet,nrow=1) +
  scale_fill_manual(values=fill_class,guide="none") +
  scale_colour_manual(values=c("New PET discordance"=white,
                               "Resolved PET discordance"=white,
                               "Unchanged classification"=ink),guide="none") +
  scale_x_discrete(drop=FALSE,limits=c("-","+"),position="bottom") +
  scale_y_discrete(drop=FALSE,limits=c("+","-")) +
  coord_equal() +
  labs(title="B",x="Switched assay result",y="C2N result") +
  fig_theme +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        strip.text=element_text(size=8,colour=ink,face="bold"),
        panel.spacing=grid::unit(2.8,"mm"),
        legend.position="none",
        axis.title=element_text(size=8),
        axis.text=element_text(size=8),
        plot.margin=margin(3,6,4,4))
range_abs <- max(abs(pair_flow$net_change_pp))
if (!is.finite(range_abs) || range_abs==0) range_abs <- 1
c_data <- pair_flow %>% mutate(
  source=factor(source,levels=rev(ASSAYS)),
  destination=factor(destination,levels=ASSAYS),
  annotate=sprintf("%+.1f",net_change_pp))
pC <- ggplot(c_data,aes(destination,source,fill=net_change_pp)) +
  geom_tile(colour="white",linewidth=.8) +
  geom_text(aes(label=annotate),size=2.95,colour=ink) +
  scale_x_discrete(drop=FALSE,labels=DISPLAY) +
  scale_y_discrete(drop=FALSE,labels=DISPLAY) +
  scale_fill_gradient2(low="#A6BCCB",mid="#F3F5F6",high="#C97B7D",
                       midpoint=0,limits=c(-range_abs,range_abs),
                       name="Target - source (pp)") +
  coord_equal() +
  labs(title="C",x="Destination assay",y="Source assay") +
  fig_theme +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        legend.position="bottom",legend.title=element_text(size=7),
        legend.text=element_text(size=7),
        axis.text.x=element_text(size=7.2)) +
  guides(fill=guide_colourbar(title.position="top",
                              barwidth=grid::unit(27,"mm"),
                              barheight=grid::unit(2.3,"mm")))
d_points <- pairs %>% filter(source==SOURCE_BD) %>%
  mutate(destination=factor(destination,levels=dest,labels=DISPLAY[dest]),
         classification_status=factor(classification_status,
           levels=c("Unchanged PET concordance","Resolved PET discordance",
                    "New PET discordance")))
if (nrow(d_points)!=3L*nrow(persons) ||
    anyDuplicated(d_points[c("RID","destination")]))
  stop("D: expected exactly one point per person per destination facet.",call.=FALSE)
write_csv(d_points %>% select(RID,destination,first_cl,last_cl,
                             first_pet,last_pet,A_first,A_last,B_last,
                             classification_status),
          "Figure_4D_Centiloid_individual_source_data.csv")
cl_range <- range(c(persons$first_cl,persons$last_cl, PET_CUTOFF),finite=TRUE)
cl_span <- diff(cl_range)
if (cl_span==0) cl_span <- 1
cl_lim <- c(cl_range[1] - .03*cl_span, cl_range[2] + .03*cl_span)
cl_breaks <- c(-100,-25,0,20,100,300,500)
cl_breaks <- cl_breaks[cl_breaks>=cl_lim[1] & cl_breaks<=cl_lim[2]]
if (!(20 %in% cl_breaks)) cl_breaks <- sort(unique(c(cl_breaks,20)))
cl_transform <- scales::pseudo_log_trans(sigma=25, base=10)
legend_cols <- c("New PET discordance"=red,
                 "Resolved PET discordance"=blue,
                 "Unchanged PET concordance"=grey)
d_points <- d_points %>% arrange(classification_status)
pD <- ggplot(d_points,aes(first_cl,last_cl)) +
  geom_abline(intercept=0,slope=1,linetype="dotted",colour=pale,linewidth=.35) +
  geom_vline(xintercept=PET_CUTOFF,linetype="dashed",colour="#8995A0",linewidth=.37) +
  geom_hline(yintercept=PET_CUTOFF,linetype="dashed",colour="#8995A0",linewidth=.37) +
  geom_point(aes(colour=classification_status,
                 alpha=classification_status),size=1.20,stroke=.10) +
  facet_wrap(~destination,nrow=1) +
  scale_colour_manual(values=legend_cols,breaks=names(legend_cols),
                      labels=c("New mismatch","Resolved mismatch","No change"),
                      name=NULL) +
  scale_alpha_manual(values=c("Unchanged PET concordance"=.27,
                              "Resolved PET discordance"=.92,
                              "New PET discordance"=.92),guide="none") +
  scale_x_continuous(trans=cl_transform,breaks=cl_breaks,limits=cl_lim,
                     guide=guide_axis(check.overlap=TRUE)) +
  scale_y_continuous(trans=cl_transform,breaks=cl_breaks,limits=cl_lim,
                     guide=guide_axis(check.overlap=TRUE)) +
  coord_equal() +
  labs(title="D",x="First PET (Centiloids)",y="Last PET (Centiloids)") +
  fig_theme +
  theme(panel.spacing=grid::unit(2.2,"mm"),
        axis.text=element_text(size=7),axis.title=element_text(size=8),
        strip.text=element_text(size=8,face="bold"),
        legend.position="bottom",legend.text=element_text(size=7),
        legend.key.size=grid::unit(2.9,"mm"),
        plot.margin=margin(3,6,3,4)) +
  guides(colour=guide_legend(nrow=1,byrow=TRUE,override.aes=list(alpha=1,size=1.6)))
fig <- (pA | pB) / (pC | pD) +
  plot_layout(widths=c(.84,1.16), heights=c(.82,1.18))
base <- file.path(out,"Figure_4_PET_anchored_v2")
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggsave(paste0(base,".pdf"), fig, width=FIG_WIDTH_MM,height=FIG_HEIGHT_MM,
       units="mm",device=pdf_device,bg="white",limitsize=FALSE)
ggsave(paste0(base,".png"), fig, width=FIG_WIDTH_MM,height=FIG_HEIGHT_MM,
       units="mm",dpi=SAVE_DPI,bg="white",limitsize=FALSE)
ggsave(paste0(base,".tiff"),fig, width=FIG_WIDTH_MM,height=FIG_HEIGHT_MM,
       units="mm",dpi=SAVE_DPI,device="tiff",compression="lzw",
       bg="white",limitsize=FALSE)
selection_audit <- tibble(
  total_core_participants=n_distinct(d$RID),
  paired_pet_90_n=nrow(persons),
  PET_0to0_n=transition_n[["PET 0->0"]],
  PET_0to1_n=transition_n[["PET 0->1"]],
  PET_1to1_n=transition_n[["PET 1->1"]],
  PET_1to0_n=transition_n[["PET 1->0"]],
  C2N_first_correct_PET_0to0_n=focus_n,
  C2N_other3_mean_delta_pp=focus_avg$mean_net_pp,
  C2N_other3_ci_low=focus_avg$ci_low,
  C2N_other3_ci_high=focus_avg$ci_high,
  bootstrap_replicates=BOOT_N)
write_csv(selection_audit,"Figure_4_selection_and_estimands_audit.csv")
caption <- paste0(
  "Figure 4. PET-anchored apparent plasma classification changes under counterfactual assay substitution. ",
  "A, first and last blood collections with paired amyloid PET scans (n = ",nrow(persons),") within +/-",PET_WINDOW,
  " days, defining PET-stable negative, PET-converting, PET-stable positive, and PET-positive-to-negative groups. ",
  "B, paired 2 x 2 last-draw Youden classifications in PET-stable negative participants whose first C2N plasma classification was negative: ",
  "continuing C2N (rows) versus substituting ALZpath, Fujirebio, or Janssen (columns). ",
  "Red cells indicate additional PET-discordant positive classifications; blue cells indicate resolved PET-discordant positivity. ",
  "Numbers in cells are participants (n = ",focus_n,"); the same participants are reused across the three destination comparisons. ",
  "Facet labels identify the destination assays; the three directed differences are displayed in panel C. ",
  "The mean destination-averaged change for C2N was ",sprintf("%+.1f",focus_avg$mean_net_pp),
  " percentage points (95% CI ",sprintf("%+.1f",focus_avg$ci_low),
  " to ",sprintf("%+.1f",focus_avg$ci_high),"). ",
  "This CI is a participant-level percentile bootstrap (",BOOT_N,
  " replicates), averaging the three destination contrasts within participant, conditional on frozen cross-fitted decisions. ",
  "C, corresponding target-minus-source positivity-rate differences (percentage points) for all 12 directed assay pairs; ",
  "eligibility is defined separately for each source assay, so inverse pairs need not cancel after conditioning. ",
  "D, actual first versus last PET Centiloid values for participants with paired PET (n = ",nrow(persons)," per destination), shown for C2N-to-other-assay comparisons. ",
  "Color indicates whether the destination changes PET-referenced classification concordance at the last draw. ",
  "Each facet contains one point per participant. Dashed lines mark Centiloid 20; the dotted line represents identity. ",
  "The axes use a monotone pseudo-logarithmic display transformation (sigma=25) to show values near and far from the threshold ",
  "while ticks indicate original Centiloid units. ",
  "All comparisons represent hypothetical assay substitution on the same last-draw specimen, not observed clinical switching; ",
  "PET-reference mismatch does not establish clinical misdiagnosis. Selection of C2N for detailed illustration is exploratory."
)
writeLines(caption,file.path(out,"Figure_4_caption_v2.txt"),useBytes=TRUE)
cat("\nFigure 4 saved to:",out,"\n")
print(selection_audit,width=Inf)
print(pair_flow %>% filter(source==SOURCE_BD) %>%
        select(destination,n,new_positive,resolved_positive,net_change_pp))
cat("\nNOTE: Panel D uses ALL paired-PET cases; B/C use PET 0->0 AND first-source-negative.\n")
cat("NOTE: Net switching effects are descriptive; platform selection is exploratory.\n")
