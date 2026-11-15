# FIGURE 1 | Patient enrollment flowchart

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


# Figure 1: R reconstruction of the enrollment flowchart from reported counts.
# Screening exclusions cannot be derived from the enrolled-patient workbook.
require_packages("svglite")
out<-output_folder("figure1");d<-load_clinical()
c<-read.csv(file.path(project_dir,"data/screening_counts.csv"));n<-setNames(c$N,c$Stage)
stopifnot(n["Screened"]-n["Age_under_18"]-n["Lost_followup_or_withdrawal"]==n["Enrolled"],n["Enrolled"]==nrow(d),n["TTV_positive"]==sum(d$TTV),n["TTV_negative"]==sum(d$TTV==0))
draw_flow<-function() {
 grid::grid.newpage()
 box<-function(x,y,w,h,label,size=12,just="centre") {
  grid::grid.rect(x,y,w,h,gp=grid::gpar(fill="white",col="#323232",lwd=1.2))
  grid::grid.text(label,x=if(just=="left")x-w/2+.016 else x,y=y,just=just,gp=grid::gpar(fontfamily="sans",fontsize=size,col="#323232",lineheight=1.18))
 }
 line<-function(x0,y0,x1,y1,arrow=FALSE)grid::grid.lines(c(x0,x1),c(y0,y1),gp=grid::gpar(col="#323232",lwd=1.1),arrow=if(arrow)grid::arrow(length=grid::unit(2.5,"mm"),type="closed")else NULL)
 line(.40,.88,.40,.43,TRUE);line(.40,.665,.54,.665,TRUE)
 line(.40,.31,.40,.25);line(.19,.25,.66,.25);line(.19,.25,.19,.12,TRUE);line(.66,.25,.66,.12,TRUE)
 box(.40,.93,.40,.09,sprintf("Total patients (n = %d)",n["Screened"]))
 box(.76,.665,.44,.30,sprintf("Excluded:\n1. Age <18 years (n = %d)\n2. Lost to follow-up or withdrawal\n   of treatment within 28 days\n   after ICU admission (n = %d)",n["Age_under_18"],n["Lost_followup_or_withdrawal"]),10.5,"left")
 box(.40,.37,.34,.12,sprintf("Enrolled patients\n(n = %d)",n["Enrolled"]))
 box(.19,.073,.34,.085,sprintf("TTV positive (n = %d)",n["TTV_positive"]))
 box(.66,.073,.36,.085,sprintf("TTV negative (n = %d)",n["TTV_negative"]))
}
prefix<-file.path(out,"figure1")
pdf(paste0(prefix,".pdf"),width=9,height=6.2,family="Helvetica",useDingbats=FALSE);draw_flow();dev.off()
png(paste0(prefix,".png"),width=9,height=6.2,units="in",res=300,type="cairo",bg="white");draw_flow();dev.off()
tiff(paste0(prefix,".tiff"),width=9,height=6.2,units="in",res=600,type="cairo",compression="lzw",bg="white");draw_flow();dev.off()
svglite::svglite(paste0(prefix,".svg"),width=9,height=6.2,bg="white");draw_flow();dev.off()
write_csv(c,file.path(out,"figure1_source_counts.csv"))
