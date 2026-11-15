# SUPPLEMENTARY FIGURE 1 | Adjusted subgroup forest plot

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


# Supplementary Figure 1: adjusted subgroup ORs from the interaction models.
# Original six-column forestplot, log OR axis, complete CIs, no overall-row label.
require_packages(c("forestploter","dplyr","jsonlite","svglite"))
suppressPackageStartupMessages({library(forestploter);library(grid);library(dplyr)})
ans<-subgroup_analysis(load_clinical());res<-ans$results
res$Subgroup<-gsub("_"," ",res$Subgroup,fixed=TRUE)
res$Estimator<-"Adjusted MLE interaction model"
stopifnot(nrow(res)==16L,res$Subgroup[1]=="")
estimated <- is.finite(res$OR)
stopifnot(sum(estimated) == 11L,
          all(is.finite(res$Lower[estimated])),
          all(is.finite(res$Upper[estimated])),
          all(res$Lower[estimated] > 0),
          all(res$Lower[estimated] < res$OR[estimated]),
          all(res$Upper[estimated] > res$OR[estimated]),
          all(res$Lower[estimated] >= 0.1),
          all(res$Upper[estimated] <= 8),
          !any(grepl("Firth", res$Estimator, ignore.case = TRUE)))
interaction_rows <- which(is.finite(res$Pint))
stopifnot(length(interaction_rows) == 5L,
          isTRUE(all.equal(p.adjust(res$Pint[interaction_rows], method = "BH"),
                           res$q_interaction_BH[interaction_rows], tolerance = 1e-12)))
for (j in interaction_rows) {
  strata <- res[c(j + 1L, j + 2L), ]
  stopifnot(sum(strata$N_noTTV) == res$N_noTTV[1],
            sum(strata$E_noTTV) == res$E_noTTV[1],
            sum(strata$N_TTV) == res$N_TTV[1],
            sum(strata$E_TTV) == res$E_TTV[1])
}
res$No_TTV_txt <- ifelse(is.na(res$N_noTTV), '', sprintf('%d/%d',res$E_noTTV,res$N_noTTV))
res$TTV_txt <- ifelse(is.na(res$N_TTV), '', sprintf('%d/%d',res$E_TTV,res$N_TTV))
res$OR_CI <- ifelse(is.na(res$OR), '', sprintf('%.2f (%.2f-%.2f)',res$OR,res$Lower,res$Upper))
res$Pint_txt <- ifelse(is.na(res$Pint), '', sprintf('%.3f',res$Pint))
out_dir <- output_folder("supplementary_figure1")
audit_dir <- file.path(out_dir,"audit")
dir.create(audit_dir,recursive=TRUE,showWarnings=FALSE)
write_csv <- function(x,fn) write.csv(x,file.path(out_dir,fn),row.names=FALSE,fileEncoding="UTF-8")
stopifnot(nrow(res)==16L, res$N_noTTV[1]==440L, res$E_noTTV[1]==231L,
          res$N_TTV[1]==42L, res$E_TTV[1]==24L)
tmp <- res |>
  dplyr::select(Subgroup, No_TTV = No_TTV_txt, TTV = TTV_txt, OR_CI, Pint_txt)
tmp$plot_col <- '        '
plot_data <- tmp |>
  dplyr::select(Subgroup, No_TTV, TTV, plot_col, OR_CI, Pint_txt)
colnames(plot_data) <- c('Subgroup', 'No TTV', 'TTV', '        ',
                         'aOR (95% CI)', 'P for interaction')
plot_data$`        ` <- paste(rep(' ', 30), collapse = '')
write_csv(plot_data, 'forest_plot_display_table.csv')

## Six-column forestploter layout retained from the supplied code.
## Reviewer changes: logarithmic axis, complete CIs, adjusted-OR header.
## Single forest plot; comparable multi-panel alignment is not applicable.
p <- forest(
  data = plot_data,
  est = res$OR,
  lower = res$Lower,
  upper = res$Upper,
  ci_column = 4,
  ref_line = 1,
  theme = forestploter::forest_theme(base_size = 12, base_family = 'sans'),
  xlim = c(0.1, 8),
  x_trans = "log",
  ticks_at = c(0.125, 0.5, 1, 2, 8),
  xlab = "Adjusted odds ratio (log scale)",
  line_height = grid::unit(7, "mm"),
  panel_space = 1.1,
  text_size = 0.8,
  colwidth = grid::unit.c(
    grid::unit(6.5, "cm"),   # Subgroup
    grid::unit(3.0, "cm"),   # No TTV
    grid::unit(3.0, "cm"),   # TTV
    grid::unit(10,  "cm"),   # Forest column
    grid::unit(5.5, "cm"),   # OR (95% CI)
    grid::unit(3.0, "cm")    # P for interaction
  ),
  point_size = 3,
  ci_lty = 1,
  ci_lwd = 1.5,
  ref_lty = 2,
  ref_lwd = 1,
  ref_col = "grey40"
)

## 6. Export by the measured grob size, not a large fixed canvas -------------
## This changes only the surrounding canvas, not the above plotting call.
## Base PDF uses explicit text sizes, avoiding Cairo's 1-point font / text-
## matrix encoding, which some PDF text-audit tools misread. Plotmath supplies
## the same >= glyph; a typographic dash becomes a hyphen only in the PDF.
p_pdf <- p
for (j in seq_along(p_pdf$grobs)) {
  gg <- p_pdf$grobs[[j]]
  if (inherits(gg, 'text') && is.character(gg$label)) {
    if (length(gg$label) == 1L && startsWith(gg$label, '\u2265')) {
      gg$label <- parse(text = paste0('phantom() >= ', substring(gg$label, 2L)))
    } else {
      gg$label <- gsub('\u2013', '-', gg$label, fixed = TRUE)
    }
    p_pdf$grobs[[j]] <- gg
  }
}
measure_pdf <- file.path(audit_dir, 'measurement_device.pdf')
grDevices::pdf(measure_pdf, width = 18, height = 12, family = 'Helvetica',
               useDingbats = FALSE)
grid::grid.newpage()
grid::grid.draw(p_pdf)
content_width <- grid::convertWidth(sum(p_pdf$widths), 'in', valueOnly = TRUE)
content_height <- grid::convertHeight(sum(p_pdf$heights), 'in', valueOnly = TRUE)
dev.off()
pad_in <- 0.06
width_in <- content_width + 2 * pad_in
height_in <- content_height + 2 * pad_in
stopifnot(is.finite(width_in), is.finite(height_in), width_in > 0, height_in > 0)

draw_tight <- function(grob = p) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = grid::unit(content_width, 'in'),
                                   height = grid::unit(content_height, 'in')))
  grid::grid.draw(grob)
  grid::popViewport()
}
prefix <- file.path(out_dir, 'supplementary_figure1')
grDevices::pdf(paste0(prefix, '.pdf'), width = width_in,
               height = height_in, family = 'Helvetica', useDingbats = FALSE, bg = 'white')
draw_tight(p_pdf); dev.off()
grDevices::png(paste0(prefix, '.png'), width = width_in, height = height_in,
              units = 'in', res = 300, type = 'cairo', bg = 'white')
draw_tight(); dev.off()
grDevices::tiff(paste0(prefix, '.tiff'), width = width_in, height = height_in,
               units = 'in', res = 600, type = 'cairo', compression = 'lzw', bg = 'white')
draw_tight(); dev.off()
if (requireNamespace('svglite', quietly = TRUE)) {
  svglite::svglite(paste0(prefix, '.svg'), width = width_in, height = height_in, bg = 'white')
  draw_tight(); dev.off()
}


# Optional Cairo PDF for applications requiring embedded TrueType fonts.
# The standard PDF above preserves explicit text sizes and Symbol glyphs.
if (identical(Sys.getenv("TTV_OPTIONAL_CAIRO_PDF", "false"), "true")) {
  grDevices::cairo_pdf(paste0(prefix, "_embedded_fonts.pdf"), width = width_in,
                       height = height_in, family = "Helvetica", bg = "white")
  draw_tight(); dev.off()
}

jsonlite::write_json(list(axis="log",limits=c(.1,8),minimum_CI=min(res$Lower,na.rm=TRUE),maximum_CI=max(res$Upper,na.rm=TRUE),complete_cases=482,all_CIs_in_frame=TRUE),file.path(audit_dir,"render_metadata.json"),pretty=TRUE,auto_unbox=TRUE)
utils::write.csv(res,file.path(out_dir,"supplementary_figure1_source_data.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
writeLines(c("Supplementary Figure 1. Subgroup analyses of the association between CSF TTV detection and the composite unfavorable outcome.","Outcome: death or coma at hospital discharge. No TTV and TTV columns show unfavorable-outcome events / patients. Analysis includes 482 complete cases.","Adjusted ORs compare TTV-positive with TTV-negative patients. Continuous age and severity are retained within strata; GCS replaces SOFA in the GCS model.","Interaction P values are unadjusted likelihood-ratio tests on identical patients. BH-FDR q values across five tests are reported in Supplementary Table 5. GCS and SOFA were assessed at ICU admission."),file.path(out_dir,"supplementary_figure1_legend.txt"))
