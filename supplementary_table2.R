# SUPPLEMENTARY TABLE 2 | CSF parameters

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


# Supplementary Table 2: eight CSF variables; missing values excluded per test.
d<-load_clinical();out<-output_folder("supplementary_table2")
r<-comparison_table(d,csf_vars,csf_labels,FALSE,TRUE);r$Unit<-csf_units
write_table(display_comparison(r),r,out,"supplementary_table2",c("Data are median [IQR]. Two-sided Wilcoxon rank-sum tests, exact=FALSE and correct=FALSE. BH-FDR across all eight CSF comparisons.","Units: chloride and glucose, mmol/L; protein, g/L; WBC and RBC, cells/microliter; MNC and PMN, %; pressure, mmH2O. Results were from the first CSF examination during hospitalization.","Full-precision CSV includes the available N in each group for each variable."))
