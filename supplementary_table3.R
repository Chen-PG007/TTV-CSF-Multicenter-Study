# SUPPLEMENTARY TABLE 3 | Imaging sensitivity analysis

# Run from RStudio (open TTV_final.Rproj first), or Rscript this_file.R.
.entry_files <- unlist(lapply(sys.frames(),function(e) if(exists("ofile",e,inherits=FALSE)) get("ofile",e) else NULL),use.names=FALSE)
.file_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
.entry <- if(length(.entry_files)) tail(.entry_files,1) else if(length(.file_arg)) sub("^--file=","",.file_arg[1]) else NULL
.candidates <- unique(c(if(!is.null(.entry)) dirname(normalizePath(.entry,mustWork=FALSE)),getwd()))
.candidates <- unique(c(.candidates,dirname(.candidates)))
.found <- .candidates[file.exists(file.path(.candidates,".ttv-project"))]
if(!length(.found)) stop("Open TTV_final.Rproj, then source this script.")
project_dir <- normalizePath(.found[1],winslash="/",mustWork=TRUE)
source(file.path(project_dir,"support/00_common.R"),local=TRUE)


# Supplementary Table 3: exact conditional imaging analysis, stratified by craniotomy.
# The exact conditional MLE is used, not the asymptotic Mantel-Haenszel estimate.
d<-load_clinical();out<-output_folder("supplementary_table3");r<-imaging_comparisons(d);arrays<-list();strata<-list()
for(i in seq_along(imaging_vars)) {
 v<-imaging_vars[i];z<-d[complete.cases(d[,c("TTV","Craniotomy_before_ICU_admission",v)]),]
 tab<-table(factor(z$TTV,levels=0:1),factor(z[[v]],levels=0:1),factor(z$Craniotomy_before_ICU_admission,levels=0:1))
 fit<-mantelhaen.test(tab,exact=TRUE,alternative="two.sided",conf.level=.95)
 r$Common_OR[i]<-unname(fit$estimate);r$Lower[i]<-fit$conf.int[1];r$Upper[i]<-fit$conf.int[2];r$Adjusted_exact_P[i]<-fit$p.value
 arrays[[v]]<-tab
 strata[[v]]<-do.call(rbind,lapply(0:1,function(k) {q<-z[z$Craniotomy_before_ICU_admission==k,];data.frame(Variable=v,Craniotomy=k,N_negative=sum(q$TTV==0),E_negative=sum(q[[v]][q$TTV==0]),N_positive=sum(q$TTV==1),E_positive=sum(q[[v]][q$TTV==1]))}))
}
r$Adjusted_BH_q<-p.adjust(r$Adjusted_exact_P,"BH")
disp<-r[,c("Modality","Finding","Overall","TTV_negative","TTV_positive")]
disp$`Craniotomy-adjusted common OR (95% CI)`<-sprintf("%.2f (%.2f-%.2f)%s",r$Common_OR,r$Lower,r$Upper,ifelse(r$Positive_events==0,"*",""))
disp$`Adjusted exact P value`<-fmt_p(r$Adjusted_exact_P);disp$`q value (BH-FDR)`<-fmt_p(r$Adjusted_BH_q)
write_table(disp,r,out,"supplementary_table3",c("Common ORs compare TTV-positive with TTV-negative patients, stratified by craniotomy before ICU admission. ORs are exact conditional maximum-likelihood estimates; 95% CIs and two-sided P values are exact conditional results from mantelhaen.test(exact=TRUE).","BH-FDR correction includes all 16 imaging comparisons (8 CT, 8 MRI). An asterisk marks zero events in the TTV-positive group. Missing findings are excluded feature by feature."))
write_csv(do.call(rbind,strata),file.path(out,"Stratified_counts.csv"));saveRDS(arrays,file.path(out,"Exact_2x2x2_tables.rds"))
