# SUPPLEMENTARY TABLE 4 | Outcome-model sensitivity analyses

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


# Supplementary Table 4: model sensitivity analyses, all ordinary logistic models.
d<-model_data(load_clinical());out<-output_folder("supplementary_table4")
fs<-list(primary_formula,update(primary_formula,.~.-Center),update(primary_formula,.~.-Steroids),update(primary_formula,.~.-SOFA_score_at_ICU_admission_+GCS_ICU))
labels<-c("Primary model","Without center fixed effects","Without steroid covariate","GCS replacing SOFA")
fits<-lapply(fs,function(f)fit_logistic(f,d[complete.cases(d[,all.vars(f)]),,drop=FALSE]))
r<-do.call(rbind,lapply(seq_along(fits),function(i) {a<-wald_coefficients(fits[[i]]);a<-a[a$term=="TTV",];data.frame(Model=labels[i],N=nobs(fits[[i]]),a)}))
disp<-data.frame(Model=r$Model,n=r$N,TTV_aOR=sprintf("%.2f",r$OR),`95% Wald CI`=sprintf("%.2f-%.2f",r$Lower,r$Upper),`P value`=fmt_p(r$P_value),check.names=FALSE)
write_table(disp,r,out,"supplementary_table4",c("Outcome: death or coma at hospital discharge (1); favorable outcome (0). ORs and 95% CIs use the Wald method; P values are two-sided Wald tests.","Primary adjustment: TTV, age per 10 years, sex, documented immunosuppression, pre-ICU craniotomy, recorded steroids, admission SOFA per point, and center fixed effects (center F reference). Alternative models are indicated in the table. Each model uses complete cases for its own variables."))
saveRDS(fits,file.path(out,"Models.rds"))
