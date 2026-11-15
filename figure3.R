# FIGURE 3 | Primary adjusted logistic model

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


# Figure 3: primary adjusted logistic model; original four-column forest scheme.
require_packages(c("forestploter","jsonlite","svglite"))
suppressPackageStartupMessages({library(forestploter);library(grid)})
out_dir<-output_folder("figure3");audit_dir<-file.path(out_dir,"audit");dir.create(audit_dir,showWarnings=FALSE)
dat<-primary_data(model_data(load_clinical()));fit_primary<-fit_logistic(primary_formula,dat)
full_coefficients<-function(fit) {
 r<-wald_coefficients(fit);r$term<-sub("^SexMale$","Sex",r$term);r$term<-sub("^Craniotomy_before_ICU_admissionYes$","Craniotomy_before_ICU_admission",r$term)
 names(r)[names(r)=="Lower"]<-"CI_low";names(r)[names(r)=="Upper"]<-"CI_high";r
}
write_csv<-function(x,filename) utils::write.csv(x,file.path(out_dir,filename),row.names=FALSE,na="",fileEncoding="UTF-8")
export_existing_main_style <- function(fit, severity_term, severity_label, stem) {
  selected <- c('TTV', 'Age_per_10_years', 'Sex', 'Immunosuppression',
                'Craniotomy_before_ICU_admission', 'Steroids', severity_term)
  labels <- c('TTV', 'Age (per 10 years)', 'Sex', 'Immunosuppression',
              'Craniotomy before ICU admission', 'Steroids', severity_label)
  coefficients <- full_coefficients(fit)
  cc <- coefficients[match(selected, coefficients$term), , drop = FALSE]
  stopifnot(!anyNA(cc$OR))
  plot_results <- data.frame(Variables = labels, OR = cc$OR,
    CI_low = cc$CI_low, CI_high = cc$CI_high, P_value = cc$P_value,
    stringsAsFactors = FALSE)
  plot_results[['aOR (95% CI)']] <- sprintf('%.2f (%.2f-%.2f)',
    plot_results$OR, plot_results$CI_low, plot_results$CI_high)
  plot_results[['P value']] <- fmt_p(plot_results$P_value)
  plot_table <- data.frame(
    Variables = plot_results$Variables,
    ' ' = rep('                              ', nrow(plot_results)),
    'aOR (95% CI)' = plot_results[['aOR (95% CI)']],
    'P value' = plot_results[['P value']], check.names = FALSE)
  zebra_fill <- rep(c('#EDF1F1', '#FFFFFF'), length.out = nrow(plot_table))
build_forest_plot <- function(font_family = "Times New Roman") {
  figure_style <- forest_style(
    base_size = 10.5,
    base_family = font_family,
    ci = gpar(col = "black", fill = "black", lwd = 1),
    ci_pch = 15,
    ref_line = gpar(col = "#777777", lty = "dashed", lwd = 1),
    xaxis = gpar(col = "black", fontsize = 9.5, fontfamily = font_family),
    xlab = gpar(col = "black", fontsize = 10.5, fontfamily = font_family),
    body = gpar(fill = zebra_fill, col = "black", fontfamily = font_family),
    header = gpar(fill = "white", col = "black", fontface = "bold",
                  fontfamily = font_family),
    fit = "none",
    core = list(
      fg_params = list(fontfamily = font_family, fontsize = 10.5),
      bg_params = list(col = NA)
    ),
    colhead = list(
      fg_params = list(fontfamily = font_family, fontsize = 10.5,
                       fontface = "bold"),
      bg_params = list(fill = "white", col = NA)
    )
  )

  p <- forest(
    data = plot_table,
    est = plot_results$OR,
    lower = plot_results$CI_low,
    upper = plot_results$CI_high,
    sizes = rep(0.40, nrow(plot_results)),
    ci_column = 2,
    ref_line = 1,
    style = figure_style
  ) |>
    set_xaxis(
      xlim = c(0.25, 6),
      ticks_at = c(0.25, 0.5, 1, 2, 4),
      ticks_digits = 2,
      x_trans = "log"
    ) |>
    set_labs(xlab = "Adjusted odds ratio (95% CI)")

  ## Keep the first header left-aligned and numeric headers centered.
  p <- edit_plot(p, part = "header", col = 1, which = "text",
                 gp = gpar(fontface = "bold", fontfamily = font_family),
                 hjust = 0, x = unit(0, "npc"))
  p <- edit_plot(p, part = "header", col = 3:4, which = "text",
                 gp = gpar(fontface = "bold", fontfamily = font_family),
                 hjust = 0.5, x = unit(0.5, "npc"))
  p
}
  p_raster <- build_forest_plot('Times New Roman')
  p_vector <- build_forest_plot('Times')
  grDevices::pdf(file.path(audit_dir, paste0(stem, '_measurement.pdf')),
    width = 18, height = 12, family = 'Times', useDingbats = FALSE)
  grid::grid.newpage(); grid::grid.draw(p_vector)
  cw <- grid::convertWidth(sum(p_vector$widths), 'in', valueOnly = TRUE)
  ch <- grid::convertHeight(sum(p_vector$heights), 'in', valueOnly = TRUE)
  dev.off()
  width <- cw + .12; height <- ch + .12
  draw <- function(g) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(width = grid::unit(cw, 'in'),
      height = grid::unit(ch, 'in')))
    grid::grid.draw(g); grid::popViewport()
  }
  pref <- file.path(out_dir, stem)
  grDevices::pdf(paste0(pref, '.pdf'), width = width, height = height,
    family = 'Times', useDingbats = FALSE, bg = 'white')
  draw(p_vector); dev.off()
  grDevices::png(paste0(pref, '.png'), width = width, height = height,
    units = 'in', res = 300, type = 'cairo', bg = 'white')
  draw(p_raster); dev.off()
  grDevices::tiff(paste0(pref, '.tiff'), width = width, height = height,
    units = 'in', res = 600, compression = 'lzw', type = 'cairo', bg = 'white')
  draw(p_raster); dev.off()
  if (requireNamespace('svglite', quietly = TRUE)) {
    svglite::svglite(paste0(pref, '.svg'), width = width, height = height,
      bg = 'white', system_fonts = list(serif = 'Times New Roman'))
    draw(p_raster); dev.off()
  }
  write_csv(plot_results, paste0(stem, '_source_data.csv'))
  jsonlite::write_json(list(figure = stem, width_in = width, height_in = height,
    padding_in = .06, raster_dpi = 600, patient_N = nobs(fit),
    panel_alignment = 'not applicable: single plot area'),
    file.path(audit_dir, paste0(stem, '_layout.json')), pretty = TRUE, auto_unbox = TRUE)
}

export_existing_main_style(fit_primary,"SOFA_score_at_ICU_admission_","SOFA score at ICU admission (per 1 point)","figure3")
write_csv(full_coefficients(fit_primary),"figure3_all_coefficients.csv");saveRDS(fit_primary,file.path(out_dir,"Primary_model.rds"))
