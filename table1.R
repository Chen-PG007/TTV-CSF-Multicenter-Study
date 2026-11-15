# TABLE 1 | Patient characteristics

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


# Table 1: 49 patient-characteristic comparisons, with one BH-FDR family.
d<-load_clinical();out<-output_folder("table1")
m<-read.csv(file.path(project_dir,"data/table1_variables.csv"),check.names=FALSE)
stopifnot(nrow(m)==49L,!any(pathogens[1:7]%in%m$source))
r<-comparison_table(d,m$source,m$variable,m$test=="Fisher exact")
r$Table_row<-m$Table_row
write_table(display_comparison(r),r,out,"table1",c("Data are n (%) or median [IQR]. Binary variables: two-sided Fisher's exact test. Continuous variables: two-sided Wilcoxon rank-sum test (normal approximation without continuity correction).","BH-FDR correction: 49 comparisons. Other viruses excludes TTV; bacteria and fungi indicate at least one detected organism in the respective category.","Units: age, years; LOS, ICU LOS and admission-to-mNGS interval, days; blood WBC, lymphocytes and neutrophils, x10^9/L; CRP, mg/L; PCT, ng/mL; SOFA, GCS and component scores, points."))
