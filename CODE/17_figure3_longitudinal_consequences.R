suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
required <- c('sb_theme_classic','sb_read_result','sb_read_core','sb_save_csv',
              'sb_require_results','sb_assert_code_current','sb_ggsave')
for (f in required) if (!exists(f, mode='function'))
  stop('Run the updated 01_config.R first; missing helper: ', f, call.=FALSE)
sb_assert_code_current()
res <- if (exists('study_b_dir')) study_b_dir else file.path(Sys.getenv('STUDY_B_ROOT'),'RESULTS')
sb_require_results('Table_S9_A3_longitudinal_pseudo_switch.csv',
                   '13_clinical_switching_and_transport.R',res)
metrics <- sb_read_result(file.path(res,'Table_S9_A3_longitudinal_pseudo_switch.csv'),check.names=FALSE)
core3 <- sb_read_core(file.path(res,'02_core_results.rds'))
rd <- core3$sb_rep_dat
assays <- c('AlzPath','C2N_ratio','Fujirebio','Janssen')
alabs <- c('ALZpath','C2N','Fujirebio','Janssen'); names(alabs) <- assays
cols <- c('RID','assay','plasma_date','sample_id','margin_youden','three_class')
if (!all(cols %in% names(rd))) stop('Core repeated data lack required columns.',call.=FALSE)
fl3 <- rd %>% mutate(assay=as.character(assay),class=as.character(three_class)) %>%
  arrange(RID,assay,plasma_date,sample_id) %>% group_by(RID,assay) %>%
  summarise(first_margin=first(margin_youden),last_margin=last(margin_youden),
            first_class=first(class),last_class=last(class),first_date=first(plasma_date),
            last_date=last(plasma_date),first_sample=first(sample_id),last_sample=last(sample_id),.groups='drop')
changes <- bind_rows(lapply(assays,function(src) bind_rows(lapply(setdiff(assays,src),function(dst) {
  s <- fl3 %>% filter(assay==src) %>% select(-assay) %>% rename_with(~paste0('s_',.x),-RID)
  t <- fl3 %>% filter(assay==dst) %>% select(-assay) %>% rename_with(~paste0('t_',.x),-RID)
  inner_join(s,t,by='RID') %>% mutate(from_assay=src,to_assay=dst)
}))))
if (anyNA(changes[c('s_first_date','s_last_date','t_first_date','t_last_date')]) ||
    !all(changes$s_first_date==changes$t_first_date) || !all(changes$s_last_date==changes$t_last_date) ||
    !all(changes$s_first_sample==changes$t_first_sample) || !all(changes$s_last_sample==changes$t_last_sample) ||
    !all(changes$s_last_date>changes$s_first_date)) stop('A3 first/last sample alignment failed.',call.=FALSE)
changes <- changes %>% transmute(RID,from_assay,to_assay,
                                 delta_same=s_last_margin-s_first_margin,delta_switch=t_last_margin-s_first_margin,
                                 delta_platform=t_last_margin-s_last_margin,
                                 same_assay_three_class_changed=s_first_class!=s_last_class,
                                 switched_three_class_changed=s_first_class!=t_last_class) %>%
  mutate(direction_reversal=delta_same*delta_switch<0,
         platform_exceeds_same_change=abs(delta_platform)>abs(delta_same),
         category_change_excess_due_to_switch=as.numeric(switched_three_class_changed)-as.numeric(same_assay_three_class_changed))
if (any(!is.finite(as.matrix(changes[c('delta_same','delta_switch','delta_platform')]))))
  stop('Non-finite A3 changes; do not silently omit observations.',call.=FALSE)
if (anyDuplicated(changes[c('RID','from_assay','to_assay')]) ||
    any(table(changes$RID)!=12L)) stop('Expected 12 unique directed comparisons per participant.',call.=FALSE)
mn <- c('direction_reversal','platform_exceeds_same_change','same_assay_three_class_changed',
        'switched_three_class_changed','category_change_excess_due_to_switch')
used <- metrics %>% filter(metric %in% mn)
for (m in mn) {
  z <- used %>% filter(metric==m)
  if (nrow(z)!=13L || anyDuplicated(z[c('from_assay','to_assay')]) ||
      any(!is.finite(as.matrix(z[c('per_100','per_100_low','per_100_high')]))) ||
      any(z$per_100_low>z$per_100 | z$per_100_high<z$per_100))
    stop('Invalid frozen estimates for ',m,call.=FALSE)
  for (j in seq_len(nrow(z))) {
    d <- if (z$from_assay[j]=='All' && z$to_assay[j]=='All') changes else
      changes %>% filter(from_assay==z$from_assay[j],to_assay==z$to_assay[j])
    est <- d %>% group_by(RID) %>% summarise(v=mean(.data[[m]]),.groups='drop') %>% summarise(v=100*mean(v))
    if (!nrow(d) || abs(est$v-z$per_100[j])>1e-6)
      stop('Reconstructed observations disagree with frozen 13: ',m,call.=FALSE)
  }
}
ov <- used %>% filter(from_assay=='All',to_assay=='All')
getm <- function(m) ov[match(m,ov$metric),]
fmt <- function(m) { z <- getm(m); sprintf('%.1f%% (95%% CI %.1f–%.1f)',z$per_100,z$per_100_low,z$per_100_high) }
ink <- '#253444'; blue <- '#254F70'; red <- '#B44145'; pale <- '#E7EBEE'
th <- sb_theme_classic() + theme(legend.position='bottom',legend.title=element_blank(),
                                 legend.text=element_text(size=7.5),plot.title=element_text(face='bold',size=12),
                                 plot.subtitle=element_text(size=7.5,colour=ink),axis.title=element_text(size=8),
                                 axis.text=element_text(size=7.5),plot.margin=margin(5,7,5,5))
lim <- max(abs(c(changes$delta_same,changes$delta_switch))) * 1.04
pA <- ggplot(changes,aes(delta_same,delta_switch)) +
  annotate('rect',xmin=-Inf,xmax=0,ymin=0,ymax=Inf,fill=red,alpha=.065) +
  annotate('rect',xmin=0,xmax=Inf,ymin=-Inf,ymax=0,fill=red,alpha=.065) +
  geom_hline(yintercept=0,colour='#ABB4BC',linewidth=.3) +
  geom_vline(xintercept=0,colour='#ABB4BC',linewidth=.3) +
  geom_abline(slope=1,intercept=0,linetype=2,colour='#74818D',linewidth=.4) +
  geom_point(aes(colour=direction_reversal),size=.6,alpha=.23) +
  scale_colour_manual(values=c('FALSE'=blue,'TRUE'=red),guide='none') +
  coord_equal(xlim=c(-lim,lim),ylim=c(-lim,lim)) +
  labs(title='A',
       x='Same-assay change (margin units)',y='Change after switching (margin units)') + th
mag <- changes %>% select(RID,from_assay,to_assay,delta_same,delta_platform) %>%
  pivot_longer(c(delta_same,delta_platform),names_to='component',values_to='change') %>%
  mutate(magnitude=abs(change),component=factor(component,levels=c('delta_same','delta_platform'),
                                                labels=c('Same-assay change','Platform shift')))
pB <- ggplot(mag,aes(magnitude,colour=component)) + stat_ecdf(geom='step',linewidth=.75,pad=FALSE) +
  scale_colour_manual(values=c('Same-assay change'=blue,'Platform shift'=red)) +
  scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1),labels=function(x)paste0(100*x,'%')) +
  scale_x_continuous(expand=expansion(mult=c(0,.025))) +
  labs(title='B',
       x='Absolute magnitude (margin units)',y='Cumulative proportion') + th
pairs <- used %>% filter(from_assay!='All',to_assay!='All',
                         metric %in% c('direction_reversal','platform_exceeds_same_change')) %>%
  mutate(from_assay=factor(from_assay,levels=rev(assays)),to_assay=factor(to_assay,levels=assays),
         endpoint=factor(metric,levels=c('direction_reversal','platform_exceeds_same_change'),
                         labels=c('Direction reversal','Platform shift >\nsame-assay change')))
pC <- ggplot(pairs,aes(to_assay,from_assay,fill=per_100)) +
  geom_tile(colour='white',linewidth=.8) + geom_text(aes(label=sprintf('%.1f',per_100)),size=2.7,colour=ink) +
  facet_wrap(~endpoint,nrow=1) + scale_x_discrete(labels=alabs,drop=FALSE) +
  scale_y_discrete(labels=alabs,drop=FALSE) +
  scale_fill_gradient(low='#F2F4F6',high='#C46565',limits=c(0,100),name='Events per 100',
                      breaks=c(0,25,50,75,100)) + coord_equal() +
  labs(title='C',x='Follow-up assay',y='Baseline assay') + th +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        strip.background=element_blank(),strip.text=element_text(size=7.5),
        axis.text.x=element_text(angle=35,hjust=1),legend.title=element_text(size=7)) +
  guides(fill=guide_colourbar(title.position='top',barwidth=grid::unit(25,'mm'),barheight=grid::unit(2,'mm')))
catdat <- ov %>% filter(metric %in% c('same_assay_three_class_changed','switched_three_class_changed')) %>%
  mutate(condition=factor(metric,levels=c('same_assay_three_class_changed','switched_three_class_changed'),
                          labels=c('Same assay','Assay switched')))
ex <- getm('category_change_excess_due_to_switch')
pD <- ggplot(catdat,aes(condition,per_100,fill=condition)) +
  geom_col(width=.48) + geom_errorbar(aes(ymin=per_100_low,ymax=per_100_high),width=.12,linewidth=.55,colour=ink) +
  geom_text(aes(y=per_100_high,label=sprintf('%.1f%%',per_100)),vjust=-.7,size=3,colour=ink) +
  scale_fill_manual(values=c('Same assay'=blue,'Assay switched'=red),guide='none') +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25),expand=expansion(mult=c(0,0))) +
  annotate('text',x=1.5,y=87,label=sprintf('Paired excess: %+.1f percentage points\n95%% CI %.1f to %.1f',
                                           ex$per_100,ex$per_100_low,ex$per_100_high),size=3,colour=red,lineheight=1.15) +
  labs(title='D',x=NULL,y='Three-class changes per 100 comparisons') + th
fig <- (pA | pB) / (pC | pD) + plot_layout(heights=c(1,1),widths=c(1.2,1))
base <- file.path(res,'Figure_3_longitudinal_switch_consequences')
pdf_device <- if (capabilities('cairo')) grDevices::cairo_pdf else grDevices::pdf
sb_ggsave(paste0(base,'.pdf'),fig,width=190,height=165,units='mm',device=pdf_device,bg='white')
sb_ggsave(paste0(base,'.png'),fig,width=190,height=165,units='mm',dpi=600,device=grDevices::png,bg='white')
sb_ggsave(paste0(base,'.tiff'),fig,width=190,height=165,units='mm',dpi=600,device=grDevices::tiff,compression='lzw',bg='white')
sb_save_csv(changes,file.path(res,'Figure_3_distribution_source_data.csv'),row.names=FALSE,na='')
sb_save_csv(used,file.path(res,'Figure_3_source_data.csv'),row.names=FALSE,na='')
writeLines(paste0('Figure 3. Distortion of apparent longitudinal change under counterfactual assay switching. ',
                  'A, same-assay versus switched change in the frozen Youden decision margin; each point is a participant and directed assay pair. ',
                  'Shaded opposite-sign quadrants identify strict direction reversals; the dashed identity line indicates equal changes. ',
                  'B, empirical cumulative distributions of absolute same-assay change and absolute follow-up platform shift. ',
                  'Platform shift equals switched change minus same-assay change. These are dependent descriptive observations, not independent tests. ',
                  'All 12 directed pairs contribute equally within every participant; neither panel omits tail observations. ',
                  'C, directed-pair reversal and platform-dominance rates per 100 comparisons; diagonal cells are inapplicable. ',
                  'D, three-class change rates with existing participant-cluster bootstrap 95% CIs. The annotated paired excess ',
                  'has its own frozen CI, not a difference of marginal interval bounds. ',
                  'All rates and CIs come from script 13, and all reconstructed point estimates are checked against it. ',
                  'Margin units use the frozen calibration and are not raw concentrations or annualized rates. ',
                  'Same-assay change is a measured change, not a direct measure of true biological progression. ',
                  'The comparisons represent same-draw counterfactual switches, not observed clinical switching.'),
           file.path(res,'Figure_3_caption.txt'))
message('Figure 3 saved; frozen point estimates matched for all 65 metric/pair rows. Participants: ',n_distinct(changes$RID))
if (interactive()) print(fig)
