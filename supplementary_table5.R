# SUPPLEMENTARY TABLE 5 | Subgroup interaction tests

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


# Supplementary Table 5: five ordinary logistic interaction LRTs, BH-FDR across five.
out<-output_folder("supplementary_table5");ans<-subgroup_analysis(load_clinical());r<-ans$interactions
disp<-data.frame(Subgroup=c("Sex","Age","Pre-ICU craniotomy","Admission GCS","Admission SOFA"),n=r$N_model,`Interaction P`=fmt_p(r$P_interaction),`BH-FDR q`=fmt_p(r$q_interaction_BH),check.names=FALSE)
write_table(disp,r,out,"supplementary_table5",c("Interaction P values are likelihood-ratio tests comparing nested models fitted to the same 482 primary complete cases. BH-FDR is applied to the five unrounded interaction P values.","Age: <65 vs >=65 years; GCS: <8 vs >=8; SOFA: <6 vs >=6. Continuous age and severity adjustment is retained. GCS replaces SOFA in the GCS interaction model. Subgroup ORs are obtained from the interaction models by changing the subgroup reference level; they are not separate unadjusted subgroup regressions."))
write_csv(ans$results,file.path(out,"Subgroup_ORs_and_counts.csv"));write_csv(ans$formulas,file.path(out,"Interaction_formulas.csv"));saveRDS(ans$models,file.path(out,"Interaction_models.rds"))
