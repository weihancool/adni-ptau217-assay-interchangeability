options(stringsAsFactors=FALSE)
root <- Sys.getenv('STUDY_B_RESULTS_DIR', unset='')
if(!nzchar(root)){
  cs<-unique(c(file.path(getwd(),'RESULTS'),getwd(),file.path(dirname(getwd()),'RESULTS'),file.path(dirname(dirname(getwd())),'RESULTS')))
  x<-cs[file.exists(file.path(cs,'02_core_results.rds'))]
  if(!length(x))stop('Set STUDY_B_RESULTS_DIR to your existing RESULTS directory')
  root<-x[1]
}
core<-readRDS(file.path(root,'02_core_results.rds'))
d<-as.data.frame(core$sb_scored)
required<-c('RID','sample_id','assay','primary_anchor','youden_positive','three_class',
  'amyloid_positive','calibrated_probability','triage_positive','confirmatory_positive')
if(!all(required%in%names(d)))stop('Missing columns: ',paste(setdiff(required,names(d)),collapse=', '))
d<-d[d$primary_anchor%in%TRUE,,drop=FALSE]
if(anyDuplicated(d[c('RID','assay')]))stop('One anchor per RID+assay required')
assays<-c('AlzPath','Fujirebio','Janssen')
if(!all(assays%in%d$assay))stop('Missing one of three prespecified assay names')
d3<-d[d$assay%in%assays,,drop=FALSE]
ids<-sort(unique(d3$RID))
if(any(table(d3$RID)!=3))stop('Incomplete three-assay anchors')
outdir<-file.path(root,'NEW_P2_SIMILAR_AUC');dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
wr<-function(a,n)utils::write.csv(a,file.path(outdir,n),row.names=FALSE,na='')
mk <- function(col, dd=d3, allow_na=FALSE){
  m<-matrix(NA,nrow=length(ids),ncol=length(assays),dimnames=list(as.character(ids),assays))
  for(j in seq_along(assays)){
    z<-dd[dd$assay==assays[j],c('RID',col)]
    m[match(as.character(z$RID),rownames(m)),j]<-as.character(z[[col]])
  }
  if(anyNA(m) && !allow_na)stop('Matrix incomplete: ',col)
  m
}
y<-matrix(as.integer(mk('youden_positive')),ncol=3,dimnames=list(as.character(ids),assays))
tri<-matrix(as.integer(mk('triage_positive')),ncol=3)
conf<-matrix(as.integer(mk('confirmatory_positive')),ncol=3)
three<-mk('three_class')
petmatrix<-mk('amyloid_positive',allow_na=TRUE)
pet<-suppressWarnings(as.integer(petmatrix[,1]))
for(j in 2:3){pp<-suppressWarnings(as.integer(petmatrix[,j]));
 if(any((is.na(pp)!=is.na(pet)) | (!is.na(pp)& !is.na(pet)&pp!=pet)))stop('PET differs across assays')}
pt<-mk('calibrated_probability');probs<-matrix(as.numeric(pt),ncol=3,dimnames=list(as.character(ids),assays))
flag_any<-function(m)as.integer(apply(m,1,function(a)length(unique(a))>1))
flag_extreme<-apply(three,1,function(a)all(c('negative','positive')%in%a))
metrics<-function(ix){
  k<-length(ix)
  ans<-c(youden_any3=mean(flag_any(y)[ix]),triage_any3=mean(flag_any(tri)[ix]),
    confirm_any3=mean(flag_any(conf)[ix]),three_any3=mean(flag_any(three)[ix]),
    direct_neg_pos3=mean(flag_extreme[ix]))
  ans
}
full<-metrics(seq_along(ids))
B<-as.integer(Sys.getenv('P2_BOOT_N',unset='2000'));if(is.na(B)||B<0)stop('Invalid P2_BOOT_N')
boots<-matrix(NA_real_,nrow=B,ncol=length(full),dimnames=list(NULL,names(full)))
if(B>0){set.seed(20261016);for(b in seq_len(B))boots[b,]<-metrics(sample(seq_along(ids),length(ids),TRUE))}
interval <- function(x){if(B<25)c(NA_real_,NA_real_)else as.numeric(quantile(x,c(.025,.975),na.rm=TRUE))}
out<-do.call(rbind,lapply(seq_along(full),function(j){ci<-interval(boots[,j]);
 data.frame(metric=names(full)[j],n_participants=length(ids),rate=full[j],
   CI_low=ci[1],CI_high=ci[2],bootstrap_attempted=B)}))
wr(out,'P2_01_three_assay_disagreement.csv')
indiv<-data.frame(RID=ids,PET=pet,ALZpath=y[,1],Fujirebio=y[,2],Janssen=y[,3],
  any_3_binary=flag_any(y),any_3_three_level=flag_any(three),direct_neg_pos=flag_extreme)
wr(indiv,'P2_02_participant_level_audit.csv')
pairs<-combn(1:3,2,simplify=FALSE)
pair_summary<-list();idx<-0
for(ii in pairs){
  a<-ii[1];b<-ii[2]; err<-as.integer(y[,a]!=y[,b]);
  for(st in c('all','PET_negative','PET_positive')){
    ix<-if(st=='all')which(rep(TRUE,length(ids))) else which(pet==if(st=='PET_negative')0L else 1L)
    if(!length(ix))next
    ci<-c(NA_real_,NA_real_)
    if(B>=25){set.seed(20261017+a*10+b+match(st,c('all','PET_negative','PET_positive')));
      v<-replicate(B,mean(sample(err[ix],length(ix),TRUE)));ci<-quantile(v,c(.025,.975))}
    idx<-idx+1;pair_summary[[idx]]<-data.frame(assay_1=assays[a],assay_2=assays[b],
      subgroup=st,n=length(ix),discordant_n=sum(err[ix]),discordant_rate=mean(err[ix]),
      CI_low=ci[1],CI_high=ci[2])
  }
}
wr(do.call(rbind,pair_summary),'P2_03_pairwise_disagreement.csv')
auc<-function(pred,label){ix<-which(!is.na(label)&is.finite(pred));label<-label[ix];pred<-pred[ix];
 if(length(unique(label))<2)return(NA_real_)
 r<-rank(pred,ties.method='average');n1<-sum(label==1);n0<-sum(label==0)
 (sum(r[label==1])-n1*(n1+1)/2)/(n1*n0)}
auc_table<-data.frame(assay=assays,n_pet=sum(!is.na(pet)),n_pet_positive=sum(pet==1,na.rm=TRUE),
 AUC=vapply(seq_along(assays),function(j)auc(probs[,j],pet),numeric(1)))
wr(auc_table,'P2_04_same_sample_pet_AUC.csv')
benchmark<-list()
for(ii in pairs){a<-ii[1];b<-ii[2]
  for(t in c(0L,1L)){
    ix<-which(pet==t)
    if(length(ix)<10)next
    pa<-mean(y[ix,a]);pb<-mean(y[ix,b]);obs<-mean(y[ix,a]!=y[ix,b])
    ind<-pa*(1-pb)+(1-pa)*pb
    benchmark[[length(benchmark)+1]]<-data.frame(assay_1=assays[a],assay_2=assays[b],
      PET_status=t,n=length(ix),Pr_positive_1=pa,Pr_positive_2=pb,
      observed_discordance=obs,conditional_independent_expected=ind,
      observed_minus_expected=obs-ind)
  }
}
wr(do.call(rbind,benchmark),'P2_05_PET_conditional_independence_reference.csv')
writeLines(c('P2: Fixed original OOF classifications; all three assays evaluated on identical anchor participants.',
 'The conditional-independence benchmark is ONLY illustrative. PET dichotomization does not remove',
 'continuous amyloid burden and shared within-person biology. Observed minus expected is not',
 'an estimator of discordance independent of diagnostic accuracy.',
 'CIs are participant bootstraps CONDITIONAL ON existing OOF thresholds.'),
 file.path(outdir,'P2_README_INTERPRETATION.txt'))
message('DONE P2 -> ',outdir)
