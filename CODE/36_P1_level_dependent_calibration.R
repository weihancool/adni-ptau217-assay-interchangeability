options(stringsAsFactors = FALSE)
if (!requireNamespace('lme4', quietly=TRUE)) stop('Install lme4 with install.packages("lme4")')
root <- Sys.getenv('STUDY_B_RESULTS_DIR', unset='')
if (!nzchar(root)) {
  candidates <- unique(c(file.path(getwd(),'RESULTS'),getwd(),
    file.path(dirname(getwd()),'RESULTS'),file.path(dirname(dirname(getwd())),'RESULTS')))
  hit <- candidates[file.exists(file.path(candidates,'02_core_results.rds'))]
  if (!length(hit)) stop('Set STUDY_B_RESULTS_DIR to the existing project RESULTS folder.')
  root <- hit[1]
}
if (!file.exists(file.path(root,'02_core_results.rds'))) stop('Missing 02_core_results.rds')
outdir <- file.path(root,'NEW_P1_CALIBRATION'); dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
writeout <- function(x,nm) utils::write.csv(x,file.path(outdir,nm),row.names=FALSE,na='')
core <- readRDS(file.path(root,'02_core_results.rds'))
d <- as.data.frame(core$sb_scored)
need <- c('RID','sample_id','assay','fold','z_value','margin_youden','plasma_date',
  'years_from_baseline','baseline_amyloid_group','phase')
if (!all(need %in% names(d))) stop('Cache missing: ',paste(setdiff(need,names(d)),collapse=', '))
if (!all(is.finite(d$z_value) & is.finite(d$margin_youden))) stop('Nonfinite score/margin')
if (anyDuplicated(d[c('sample_id','assay')])) stop('Duplicate sample_id + assay')
assays <- unique(as.character(d$assay)); if(length(assays)!=4L)stop('Expected exactly four assays')
if (!all(table(d$sample_id)==4L))stop('Each sample must have all four assays')
if (any(tapply(as.character(d$fold),d$RID,function(v)length(unique(v)))!=1L))
  stop('CV leakage: RID occurs in multiple original folds')
d$assay<-factor(as.character(d$assay),levels=assays)
d$RID <- as.character(d$RID); d$fold <- as.character(d$fold)
d$sample_id <- as.character(d$sample_id)
d$phase<-as.character(d$phase)
d$baseline_amyloid_group<-as.character(d$baseline_amyloid_group)
zmean <- ave(d$z_value,d$sample_id,FUN=mean)
d$proxy_all_four<-zmean
d$proxy_leave_own_out<-(4*zmean-d$z_value)/3
calibrate <- function(data,proxy_col,method) {
  result <- rep(NA_real_,nrow(data))
  diagrows<-list(); folds<-sort(unique(data$fold))
  for (f in folds) {
    tr<-data[data$fold!=f,,drop=FALSE]; te<-data[data$fold==f,,drop=FALSE]
    tr$proxy <- tr[[proxy_col]]; te$proxy<-te[[proxy_col]]
    tr$assay<-factor(tr$assay,levels=assays); te$assay<-factor(te$assay,levels=assays)
    formula <- if(method=='linear') margin_youden ~ assay*proxy else
      margin_youden ~ assay*splines::ns(proxy,df=3)
    fit<-stats::lm(formula,data=tr)
    if(anyNA(coef(fit))) stop('Rank-deficient calibration in fold ',f,' method ',method)
    n<-nrow(te)
    pred_data<-data.frame(proxy=rep(te$proxy,times=length(assays)),
      assay=factor(rep(assays,each=n),levels=assays))
    pp<-matrix(stats::predict(fit,newdata=pred_data),nrow=n,ncol=length(assays))
    x0<-stats::median(tr$proxy)
    p0<-as.numeric(stats::predict(fit,newdata=data.frame(
      proxy=rep(x0,length(assays)),assay=factor(assays,levels=assays))))
    aidx<-match(as.character(te$assay),assays)
    correction <- (pp[cbind(seq_len(n),aidx)]-rowMeans(pp)) -
      (p0[aidx]-mean(p0))
    result[which(data$fold==f)]<-te$margin_youden-correction
    diagrows[[length(diagrows)+1L]]<-data.frame(
      proxy=proxy_col,method=method,fold=f,training_rids=length(unique(tr$RID)),
      test_rids=length(unique(te$RID)),training_rows=nrow(tr),test_rows=n,
      corrected_sd=stats::sd(correction),max_abs_correction=max(abs(correction)))
  }
  if(any(!is.finite(result))) stop('Incomplete out-of-fold predictions')
  list(value=result,audit=do.call(rbind,diagrows))
}
message('P1: learning 4 out-of-fold shape corrections ...')
specs<-list(raw=d$margin_youden); audits<-list()
for (proxy in c('proxy_all_four','proxy_leave_own_out')) {
  for (method in c('linear','spline')) {
    nm<-paste(method,if(proxy=='proxy_all_four')'all4' else 'leave1out',sep='_')
    tmp<-calibrate(d,proxy,method)
    specs[[nm]]<-tmp$value;audits[[nm]]<-tmp$audit
  }
}
writeout(do.call(rbind,audits),'P1_01_crossfit_curve_audit.csv')
prep <- function(data,y) {
  data$y<-as.numeric(y);data$RID<-factor(data$RID)
  data$assay<-factor(data$assay,levels=assays)
  data$occ_id<-factor(data$sample_id)
  data$baseline_amyloid_group<-droplevels(factor(data$baseline_amyloid_group))
  data$phase<-droplevels(factor(data$phase))
  data
}
fitvar <- function(data,y,tag) {
  dd<-prep(data,y)
  fix<-c(if(nlevels(dd$baseline_amyloid_group)>1)
    'years_from_baseline*baseline_amyloid_group' else 'years_from_baseline',
    if(nlevels(dd$phase)>1) 'phase','assay')
  fo<-stats::as.formula(paste('y ~',paste(c(fix,
    '(1|RID)','(1|occ_id)','(1|RID:assay)'),collapse=' + ')))
  m<-lme4::lmer(fo,data=dd,REML=TRUE,na.action=na.fail,
    control=lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=2e5)))
  v<-as.data.frame(lme4::VarCorr(m)); getv<-function(g){t<-v$vcov[v$grp==g];if(length(t))t[1] else NA_real_}
  vals<-c(person=getv('RID'),shared_occasion=getv('occ_id'),
    person_by_assay=getv('RID:assay'),residual=getv('Residual'))
  total<-sum(vals)
  tab<-data.frame(model=tag,component=names(vals),variance=as.numeric(vals),
    share=as.numeric(vals)/total,n_participants=length(unique(dd$RID)),
    n_occasions=length(unique(dd$occ_id)),n_rows=nrow(dd),
    singular=lme4::isSingular(m),converged=is.null(m@optinfo$conv$lme4$messages))
  summ<-data.frame(model=tag,PxA=vals['person_by_assay'],
    PxA_share=vals['person_by_assay']/total,
    PxA_over_occasion=vals['person_by_assay']/vals['shared_occasion'],
    PxA_over_occasion_residual=vals['person_by_assay']/(vals['shared_occasion']+vals['residual']),
    singular=lme4::isSingular(m),n_rid=length(unique(dd$RID)))
  list(component=tab,summary=summ)
}
persistence <- function(data,y,tag){
  dev <- y-ave(y,data$sample_id,FUN=mean)
  z<-data.frame(RID=data$RID,assay=factor(data$assay,levels=assays),
    sample_id=data$sample_id,date=as.Date(data$plasma_date),
    years=data$years_from_baseline,phase=factor(data$phase),dev=dev)
  z<-z[order(z$RID,z$assay,z$date,z$sample_id),]
  fixed<-if(nlevels(z$phase)>1)dev~assay*years+phase else dev~assay*years
  z$res<-stats::residuals(stats::lm(fixed,data=z))
  key<-interaction(z$RID,z$assay,drop=TRUE)
  byid<-split(seq_len(nrow(z)),key)
  ends<-do.call(rbind,lapply(byid,function(ix){
    if(length(unique(z$sample_id[ix]))<2)return(NULL)
    c(first=z$res[ix[1]],last=z$res[ix[length(ix)]])
  }))
  r<-if(!is.null(ends)&&nrow(ends)>2) stats::cor(ends[,1],ends[,2],method='spearman') else NA_real_
  z$RID_assay<-interaction(z$RID,z$assay,drop=TRUE)
  fixed_icc <- if(nlevels(z$phase)>1) dev~assay*years+phase+(1|RID_assay) else dev~assay*years+(1|RID_assay)
  im<-lme4::lmer(fixed_icc,data=z,REML=TRUE,
    control=lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=2e5)))
  va<-as.data.frame(lme4::VarCorr(im))
  u<-va$vcov[va$grp=='RID_assay'][1];v<-va$vcov[va$grp=='Residual'][1]
  data.frame(model=tag,spearman_first_last=r,ICC=u/(u+v),
    n_first_last_pairs=if(is.null(ends))0 else nrow(ends),
    n_participants=length(unique(z$RID)))
}
message('P1: fitting original-comparable lmer for raw and 4 adjusted scales ...')
rawcomponents<-list(); sums<-list(); persists<-list()
for (nm in names(specs)) {
  x<-fitvar(d,specs[[nm]],nm)
  rawcomponents[[nm]]<-x$component;sums[[nm]]<-x$summary
  persists[[nm]]<-persistence(d,specs[[nm]],nm)
  message('  ',nm,': PxA=',round(x$summary$PxA,4),' ICC=',round(persists[[nm]]$ICC,3))
}
tab<-do.call(rbind,rawcomponents); ss<-do.call(rbind,sums)
pe<-do.call(rbind,persists)
raw_pxa<-ss$PxA[ss$model=='raw'][1]
ss$PxA_absolute_reduction_percent<-100*(1-ss$PxA/raw_pxa)
writeout(tab,'P1_02_variance_components.csv')
writeout(ss,'P1_03_headline_comparison.csv')
writeout(pe,'P1_04_persistence_comparison.csv')
writeout(data.frame(RID=d$RID,sample_id=d$sample_id,assay=as.character(d$assay),
  fold=d$fold,margin_original=d$margin_youden,
  linear_all4=specs$linear_all4,spline_all4=specs$spline_all4,
  linear_leave1out=specs$linear_leave1out,
  spline_leave1out=specs$spline_leave1out),
  'P1_05_out_of_fold_corrected_margins.csv')
B<-as.integer(Sys.getenv('P1_BOOT_N',unset='0'))
if(is.na(B)||B<0)stop('P1_BOOT_N must be nonnegative integer')
if(B>0){
  set.seed(20261015)
  ids<-unique(d$RID); groups<-split(seq_len(nrow(d)),d$RID)
  draws<-list(); cnt<-0L
  for(b in seq_len(B)){
    picks<-sample(ids,length(ids),replace=TRUE)
    ix<-unlist(groups[picks],use.names=FALSE)
    nn<-rep(seq_along(picks),times=vapply(groups[picks],length,integer(1)))
    bd<-d[ix,,drop=FALSE]
    bd$RID<-paste0('rep',nn)
    bd$sample_id<-paste0(bd$RID,'_',bd$sample_id)
    for(nm in names(specs)) {
      z<-tryCatch({
        v<-fitvar(bd,specs[[nm]][ix],nm)$summary
        data.frame(draw=b,model=nm,PxA=v$PxA,PxA_share=v$PxA_share,
          PxA_over_occasion=v$PxA_over_occasion)
      },error=function(e){
        data.frame(draw=b,model=nm,PxA=NA_real_,PxA_share=NA_real_,
          PxA_over_occasion=NA_real_)
      })
      cnt<-cnt+1L;draws[[cnt]]<-z
    }
    if(b%%max(1L,B%/%10L)==0L)message('P1 conditional bootstrap ',b,'/',B)
  }
  boot<-do.call(rbind,draws);writeout(boot,'P1_06_conditional_bootstrap_draws.csv')
  cis<-do.call(rbind,lapply(split(boot,boot$model),function(w){
    data.frame(model=w$model[1],successful=sum(is.finite(w$PxA)),attempted=nrow(w),
      PxA_ci_low=if(sum(is.finite(w$PxA))>10)quantile(w$PxA,.025,na.rm=TRUE) else NA_real_,
      PxA_ci_high=if(sum(is.finite(w$PxA))>10)quantile(w$PxA,.975,na.rm=TRUE) else NA_real_)
  }))
  writeout(cis,'P1_07_conditional_bootstrap_CI.csv')
}
writeLines(c(
 'P1 INTERPRETATION NOTES',
 '1. Variance before and after calibration is estimated on the same sample and margin scale.',
 '2. proxy_all_four uses the mean original cross-fitted z-value; the proxy includes the outcome assay itself.',
 '3. proxy_leave_own_out averages other 3 assays and is the secondary circularity check.',
 '4. Corrections are learned on other participants (original CV folds).',
 '5. Corrections remove only assay-specific slope/shape differences relative to the mean curve,',
 '   at a fixed proxy reference level. A common level response is retained.',
 '6. ICC and first-to-last correlation are sensitivity metrics, not exact refits of Table S9.',
 '7. Bootstrap CI, if enabled, conditions on frozen cross-fitted calibration functions.',
 '8. Large reduction does NOT prove clinically validated calibration or remove binary discordance.',
 '9. If results depend strongly on all-four vs leave-own-out proxy, treat this as model ambiguity.'),
 file.path(outdir,'P1_README_INTERPRETATION.txt'))
message('DONE P1 -> ',outdir)
