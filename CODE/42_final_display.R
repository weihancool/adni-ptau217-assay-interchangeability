args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else Sys.getenv("STUDY_B_ROOT", "")
if (!nzchar(root)) {
  own <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
  root <- if (length(own)) dirname(dirname(normalizePath(own))) else getwd()
}
res <- file.path(root, "RESULTS")
out <- file.path(res, "FINAL_DISPLAY")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
read_frozen <- function(name) read.csv(file.path(res,name), stringsAsFactors=FALSE,
                                      check.names=FALSE, fileEncoding="UTF-8-BOM")
ink <- "#1E2A36"; navy <- "#1F4E79"; red <- "#B23A48"
axiscol <- "#7C8A96"; band <- "#F5F7F9"; gridcol <- "#E9EDF0"
assays <- c("AlzPath", "C2N_ratio", "Fujirebio", "Janssen")
assay_names <- c("ALZpath", "C2N", "Fujirebio", "Janssen")
components <- c("person", "person_by_assay", "occasion_within_person", "residual")
component_names <- c("Participant", "Participant \u00d7 assay", "Shared collection", "Residual")
component_cols <- c(navy, red, "#8DAAC4", "#DCE2E7")
margin <- subset(read_frozen("Table_S7_A1_switching_risk_by_source_margin_CI.csv"), metric=="binary_flip")
bench <- read_frozen("Table_S7_A1_residual_variance_benchmark.csv")
bench <- bench[match(margin$source_margin_bin, bench$source_margin_bin),]
vc <- subset(read_frozen("Figure_S5_variance_components_main.csv"), audit=="crossfit_main" & component %in% components)
deviations <- read_frozen("Figure_S6_deviation_first_vs_last.csv")
rho <- read_frozen("Figure_S6_deviation_persistence_by_assay.csv")
transitions <- read_frozen("Figure_S6_deviation_sign_transition.csv")
ratios <- read_frozen("Figure_S5_variance_ratio_ci.csv")
draws <- subset(read_frozen("Figure_S5_variance_ratio_bootstrap_draws.csv"), success & is.finite(value))
modes <- read_frozen("Figure_S5_person_by_assay_conditional_modes.csv")
stopifnot(nrow(margin)==9L, !anyNA(bench$source_margin_bin),
          max(abs(margin$per_100-bench$observed_per100))<1e-8,
          nrow(deviations)==1572L, length(unique(deviations$RID))==393L,
          nrow(modes)==1572L, all(rho$n_participants==393L))
for (s in unique(vc$scale)) stopifnot(abs(sum(vc$proportion[vc$scale==s])-1)<1e-8)
for (a in assays) {
  d <- subset(deviations, assay==a)
  r <- subset(rho, assay==a)
  z <- subset(transitions, assay==a)
  stopifnot(nrow(d)==393L, abs(cor(d$first_deviation,d$last_deviation,method="spearman")-r$deviation_spearman)<1e-10,
            sum(z$n)==393L,
            abs(mean((d$first_deviation<0)==(d$last_deviation<0))-z$sign_retention[1])<1e-10)
}
setup <- function(fig,mar,new=TRUE) {
  par(fig=fig, mar=mar, new=new, family="Arial", ps=8,
      cex=1, cex.axis=.95, cex.lab=1, col=ink, col.axis=ink, col.lab=ink,
      mgp=c(2,.45,0), tcl=-.22, lend="butt", xaxs="i", yaxs="i")
}
letter <- function(x, y, label) {
  par(fig=c(0,1,0,1),mar=rep(0,4),new=TRUE)
  plot.new(); plot.window(c(0,1),c(0,1)); text(x,y,label,adj=c(0,1),font=2,cex=1.5)
}
export <- function(name,w,h,draw) {
  for (ext in c("pdf","tiff","png")) {
    path <- file.path(out,paste0(name,".",ext))
    if (ext=="pdf") cairo_pdf(path,width=w/25.4,height=h/25.4,family="Arial",pointsize=8)
    if (ext=="tiff") tiff(path,width=w,height=h,units="mm",res=600,compression="lzw",pointsize=8,type="cairo")
    if (ext=="png") png(path,width=w,height=h,units="mm",res=300,pointsize=8,type="cairo")
    draw(); dev.off()
  }
}
draw_main <- function() {
  setup(c(0,.54,.645,1),c(5.4,6.4,2.6,.2),FALSE)
  plot(NA,xlim=c(0,75),ylim=c(.5,10),axes=FALSE,xlab="",ylab="")
  y <- 9:1
  for (j in seq(1,9,2)) rect(0,j-.5,75,j+.5,col=band,border=NA)
  abline(v=seq(0,50,10),col=gridcol,lwd=.6)
  segments(margin$per_100_low,y+.10,margin$per_100_high,y+.10,col=navy,lwd=1.25)
  points(margin$per_100,y+.10,pch=21,bg="white",col=navy,cex=.9,lwd=1)
  points(bench$residual_null_per100,y-.15,pch=18,col=red,cex=.9)
  text(margin$per_100_high+1,y+.10,sprintf("%.1f",margin$per_100),adj=0,cex=.9)
  labs <- c("0-<0.25","0.25-<0.50","0.50-<0.75","0.75-<1.00","1.00-<1.25","1.25-<1.50","1.50-<1.75","1.75-<2.00","\u22652.00")
  axis(1,at=seq(0,50,10),col=axiscol,col.ticks=axiscol)
  axis(2,at=y,labels=labs,las=1,tick=FALSE,cex.axis=.93)
  mtext("Absolute source-assay Youden margin (SD units)",side=2,line=4.2,cex=.95)
  mtext("Binary reclassifications per 100 comparisons",side=1,line=1.8,cex=.95)
  text(62,9.7,"Contributing\nparticipants",cex=.9,adj=.5,xpd=NA)
  text(62,y,margin$participants,cex=.95)
  legend("bottom",inset=c(-.13,-.40),xpd=NA,bty="n",ncol=1,cex=.95,
         legend=c("Observed","Residual-variance benchmark"),pch=c(21,18),
         pt.bg=c("white",red),col=c(navy,red),y.intersp=1.3)
  setup(c(.54,1,.645,1),c(5.4,5.6,2.6,.65))
  plot(NA,xlim=c(0,100),ylim=c(.4,4.8),axes=FALSE,xlab="",ylab="")
  scales <- c("margin_triage","margin_confirmatory","margin_youden","calibrated_probability")
  scale_labs <- c("Triage\nmargin","Confirmatory\nmargin","Youden\nmargin","Calibrated\nprobability")
  rect(0,1.52,100,2.48,col=band,border=NA)
  for (i in 1:4) {
    z <- vc[vc$scale==scales[i],]; z <- z[match(components,z$component),]
    shares <- 100*z$proportion; left <- c(0,head(cumsum(shares),-1)); yy <- 5-i
    rect(left,yy-.27,left+shares,yy+.27,col=component_cols,border="white",lwd=.5)
    for (k in 1:4) if (i==3 || shares[k]>=15) {
      label_y <- if(i==3 && k==3) yy-.45 else yy
      if(i==3 && k==3) segments(left[k]+shares[k]/2,yy-.27,left[k]+shares[k]/2,yy-.34,col=axiscol,lwd=.6)
      text(left[k]+shares[k]/2,label_y,sprintf("%.1f%%",shares[k]),cex=.9,
           col=if(k<=2) "white" else ink,font=if(i==3) 2 else 1)
    }
  }
  axis(1,at=c(0,25,50,75,100),col=axiscol,col.ticks=axiscol)
  axis(2,at=4:1,labels=scale_labs,las=1,tick=FALSE,cex.axis=.95)
  mtext("Share of random variance (%)",side=1,line=1.8,cex=.98)
  legend("bottom",inset=c(-.22,-.40),xpd=NA,bty="n",ncol=2,cex=.9,
         legend=component_names,fill=component_cols,border=NA,x.intersp=.5,y.intersp=1.3)
  lim <- range(c(deviations$first_deviation,deviations$last_deviation),finite=TRUE)
  lim <- lim+c(-1,1)*diff(lim)*.06
  frames <- list(c(.07,.52,.33,.625),c(.54,.99,.33,.625),c(.07,.52,.035,.33),c(.54,.99,.035,.33))
  for (i in 1:4) {
    setup(frames[[i]],c(3.1,3.6,2.8,.7))
    d <- subset(deviations,assay==assays[i]); r <- subset(rho,assay==assays[i])
    plot(d$first_deviation,d$last_deviation,xlim=lim,ylim=lim,asp=1,
         pch=16,cex=.48,col=adjustcolor(navy,alpha.f=.40),axes=FALSE,xlab="",ylab="")
    abline(0,1,col="#AAB2B7",lty=2,lwd=.8)
    abline(lm(last_deviation~first_deviation,data=d),col=red,lwd=1)
    axis(1,col=axiscol,col.ticks=axiscol); axis(2,las=1,col=axiscol,col.ticks=axiscol)
    mtext(assay_names[i],side=3,line=1.6,font=2,cex=1.1)
    mtext(sprintf("\u03c1 = %.3f (95%% CI %.3f-%.3f)",r$deviation_spearman,r$deviation_spearman_low,r$deviation_spearman_high),side=3,line=.35,cex=.98)
    if(i>2) mtext("First-collection adjusted deviation",side=1,line=1.8,cex=.98)
    if(i%%2==1) mtext("Last-collection adjusted deviation",side=2,line=2.25,cex=.98)
  }
  letter(.018,.987,"A");letter(.545,.987,"B");letter(.018,.632,"C")
}
draw_s3 <- function() {
  setup(c(0,1,.57,1),c(3.7,14,2.4,2),FALSE)
  plot(NA,xlim=c(0,max(draws$value,ratios$high)*1.04),ylim=c(.5,2.55),axes=FALSE,xlab="",ylab="")
  abline(v=1,col=axiscol,lty=2,lwd=.8)
  ratio_cols <- c(red,navy)
  for(i in 1:2) {
    z<-ratios[i,]; y<-3-i; values<-draws$value[draws$statistic==z$statistic]
    den<-density(values,n=512,cut=0)
    polygon(c(den$x,rev(den$x)),c(rep(y,length(den$x)),rev(y+.34*den$y/max(den$y))),
            border=NA,col=adjustcolor(ratio_cols[i],alpha.f=.3))
    segments(z$low,y-.07,z$high,y-.07,col=ratio_cols[i],lwd=1.1)
    q<-quantile(values,c(.25,.75));segments(q[1],y-.07,q[2],y-.07,col=ratio_cols[i],lwd=3)
    points(z$estimate,y-.07,pch=21,bg="white",col=ratio_cols[i],cex=1.1,lwd=1)
    text(par("usr")[2]*.99,y+.30,sprintf("%.2f (95%% CI %.2f-%.2f)",z$estimate,z$low,z$high),adj=1,cex=1)
  }
  axis(1,col=axiscol,col.ticks=axiscol)
  axis(2,at=2:1,labels=c("Participant \u00d7 assay /\nshared collection","Participant \u00d7 assay /\n(shared collection + residual)"),las=1,tick=FALSE)
  mtext("Variance ratio (Youden margin)",side=1,line=2,cex=1.05)
  setup(c(0,1,0,.54),c(5.5,5,2.5,2))
  cols<-setNames(c(navy,red,"#7FA9BC","#9A7FA8"),assays)
  plot(NA,xlim=range(modes$rank),ylim=range(modes$low,modes$high),axes=FALSE,xlab="",ylab="")
  abline(h=0,col=axiscol,lwd=.7)
  segments(modes$rank,modes$low,modes$rank,modes$high,col=adjustcolor(cols[modes$assay],alpha.f=.25),lwd=.4)
  points(modes$rank,modes$conditional_mode,col=cols[modes$assay],pch=16,cex=.22)
  axis(1,col=axiscol,col.ticks=axiscol);axis(2,las=1,col=axiscol,col.ticks=axiscol)
  mtext("Participant \u00d7 assay effects, ranked",side=1,line=2)
  mtext("Conditional mode (SD units)",side=2,line=3)
  mtext("1,572 participant-assay combinations from 393 participants",side=3,line=.6,adj=0,cex=1)
  legend("bottom",inset=c(0,-.29),xpd=NA,bty="n",horiz=TRUE,legend=assay_names,col=cols,pch=16,cex=1)
  letter(.018,.985,"A");letter(.018,.535,"B")
}
draw_s6 <- function() {
  frames<-list(c(0,.5,.5,1),c(.5,1,.5,1),c(0,.5,0,.5),c(.5,1,0,.5))
  signs<-c("Below zero","At or above zero")
  for(i in 1:4) {
    setup(frames[[i]],c(4.7,6,4.7,1.4),i!=1)
    z<-subset(transitions,assay==assays[i])
    plot(NA,xlim=c(.5,2.5),ylim=c(.5,2.5),axes=FALSE,xlab="",ylab="")
    for(j in 1:nrow(z)) {
      x<-match(z$last_sign[j],signs);y<-3-match(z$first_sign[j],signs)
      color<-if(z$first_sign[j]==z$last_sign[j]) navy else red
      p<-z$conditional_proportion[j]
      rect(x-.48,y-.48,x+.48,y+.48,col=adjustcolor(color,alpha.f=.18+.82*p),border="white",lwd=2)
      text(x,y,sprintf("%.0f%%\n(n = %d)",100*p,z$n[j]),cex=1.1,col=if(p>=.55) "white" else ink)
    }
    axis(1,at=1:2,labels=c("Below","At/above"),tick=FALSE)
    axis(2,at=2:1,labels=c("Below","At/above"),tick=FALSE,las=1)
    mtext("Last-collection deviation sign",side=1,line=2,cex=1)
    mtext("First-collection deviation sign",side=2,line=4.5,cex=1)
    mtext(assay_names[i],side=3,line=2.4,font=2,cex=1.2)
    mtext(sprintf("%.1f%% retained",100*z$sign_retention[1]),side=3,line=.9,cex=1.05)
  }
}
export("Figure_2_variance_and_persistence",190,230,draw_main)
export("Figure_S3_variance_ratios_and_conditional_modes",190,150,draw_s3)
export("Figure_S6_deviation_sign_retention",190,150,draw_s6)
old<-c("Figure_S1_flowchart","Figure_S2_gray_zone_patterns_top15","Figure_S3_threshold_distributions","Figure_S4_sensitivity_ratio","Figure_S7_remeasurement_vs_switching","Figure_S8_A3_small_change_sensitivity")
new<-c("Figure_S1_cohort_selection","Figure_S2_classification_patterns","Figure_S4_distributions_and_thresholds","Figure_S5_PET_definition_sensitivity","Figure_S7_three_class_changes","Figure_S8_longitudinal_sensitivity")
for(i in seq_along(old)) for(ext in c("pdf","tiff","png")) {
  stopifnot(file.copy(file.path(res,paste0(old[i],".",ext)),file.path(out,paste0(new[i],".",ext)),overwrite=TRUE))
}
write.csv(cbind(margin, residual_benchmark_per100=bench$residual_null_per100),file.path(out,"Figure_2A_source_data.csv"),row.names=FALSE)
write.csv(vc,file.path(out,"Figure_2B_source_data.csv"),row.names=FALSE)
write.csv(deviations,file.path(out,"Figure_2C_source_data.csv"),row.names=FALSE)
write.csv(rho,file.path(out,"Figure_2C_correlations.csv"),row.names=FALSE)
write.csv(data.frame(old_figure=c(1,2,3,4,5,6,7,8),final_figure=c(1,4,3,6,2,5,7,8)),file.path(out,"Supplementary_figure_number_map.csv"),row.names=FALSE)
write.csv(data.frame(old_table=1:14,final_table=c(7,8,10,11,12,13,1,2,3,4,9,5,6,14)),file.path(out,"Supplementary_table_number_map.csv"),row.names=FALSE)
inputs<-c("Table_S7_A1_switching_risk_by_source_margin_CI.csv","Table_S7_A1_residual_variance_benchmark.csv","Figure_S5_variance_components_main.csv","Figure_S6_deviation_first_vs_last.csv","Figure_S6_deviation_persistence_by_assay.csv","Figure_S6_deviation_sign_transition.csv","Figure_S5_variance_ratio_ci.csv","Figure_S5_variance_ratio_bootstrap_draws.csv","Figure_S5_person_by_assay_conditional_modes.csv")
write.csv(data.frame(file=inputs,md5=unname(tools::md5sum(file.path(res,inputs)))),file.path(out,"Frozen_display_input_md5.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"Display_session_info.txt"))
submission_figures <- file.path(root,"Manuscript","final","Figures")
dir.create(submission_figures,recursive=TRUE,showWarnings=FALSE)
for(ext in c("pdf","tiff","png")) {
  stopifnot(file.copy(file.path(out,paste0("Figure_2_variance_and_persistence.",ext)),
    file.path(submission_figures,paste0("Figure_2_switching_consequences_v2.",ext)),overwrite=TRUE))
}
for(path in list.files(out,pattern="^Figure_S[1-8]_.*[.](pdf|tiff|png)$",full.names=TRUE)) {
  stopifnot(file.copy(path,file.path(submission_figures,basename(path)),overwrite=TRUE))
}
message("Final submission figures exported to ",out,". All frozen point estimates and CIs preserved.")
