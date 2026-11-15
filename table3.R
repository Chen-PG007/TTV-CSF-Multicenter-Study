# TABLE 3 | CT and MRI findings

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


# Table 3: 8 CT and 8 MRI comparisons, two-sided Fisher tests, BH-FDR over 16.
d<-load_clinical();out<-output_folder("table3");r<-imaging_comparisons(d)
disp<-r[,c("Modality","Finding","Overall","TTV_negative","TTV_positive")];disp$`P value`<-fmt_p(r$P_value);disp$`q value (BH-FDR)`<-fmt_p(r$q_value)
write_table(disp,r,out,"table3",c("Data are n/N (%); N is the number with available data for the specific feature. Missing or unevaluable findings are excluded, not counted as normal.","P values: two-sided Fisher's exact tests. BH-FDR: all 16 CT and MRI comparisons. CT meningeal abnormalities and MRI meningeal enhancement are separate variables."))
