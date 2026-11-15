# SUPPLEMENTARY TABLE 1 | TTV detection by center

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


# Supplementary Table 1: detection rates at centers A-F.
d<-load_clinical();out<-output_folder("supplementary_table1")
r<-do.call(rbind,lapply(LETTERS[1:6],function(a) {x<-d[d$Center==a,];data.frame(Center=a,Total=nrow(x),TTV_positive=sum(x$TTV),Detection_rate=100*mean(x$TTV))}))
disp<-r;disp$Detection_rate<-sprintf("%.1f",r$Detection_rate)
write_table(disp,r,out,"supplementary_table1","Centers A-F are anonymized. Detection rate (%) = TTV-positive patients / all patients at each center x 100.")
