# TABLE 2 | Co-detected pathogens

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


# Table 2: 16 specific co-detected pathogens, one BH-FDR family.
d<-load_clinical();out<-output_folder("table2")
r<-comparison_table(d,pathogens,pathogen_labels,TRUE,TRUE)
r$Category<-c(rep("Viruses",7),rep("Bacteria",7),rep("Fungi",2))
write_table(display_comparison(r,FALSE),r,out,"table2",c("Data are n (%). P values are from two-sided Fisher's exact tests.","BH-FDR correction was applied across all 16 pathogen comparisons. TTV was excluded from the other-virus comparisons."))
