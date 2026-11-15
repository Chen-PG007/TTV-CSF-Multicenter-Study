# FIGURE 2 | CSF raincloud plots

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


# Figure 2: original half-violin + points + boxplot scheme, 2 rows x 4 columns.
# Protein, WBC and RBC use log10(x+1) for display only. Tests use original data.
# Stars show unadjusted P; the eight BH-FDR q values are in Supplementary Table 2.
require_packages(c("ggplot2","gghalves","e1071","dplyr","patchwork","gridExtra","jsonlite","svglite"))
suppressPackageStartupMessages({library(ggplot2);library(gghalves);library(dplyr);library(patchwork)})
d<-load_clinical();out<-output_folder("figure2")
data<-d[,c("TTV",csf_vars)];names(data)<-c("TTV",csf_labels)
data$TTV<-factor(data$TTV,levels=0:1,labels=c("TTV-","TTV+"));group_var<-"TTV"
set.seed(1963580)
make_violin_plot <- function(v, dat) {
  x_raw <- dat[[v]]
  stopifnot(all(is.na(x_raw) | is.finite(x_raw)), all(is.na(x_raw) | x_raw >= 0))
  

  if (all(is.na(x_raw)) || length(unique(na.omit(x_raw))) <= 1) {
    message("Skipping ", v, ": all values are missing or constant.")
    return(NULL)
  }
  

  sk <- e1071::skewness(x_raw, na.rm = TRUE)
  use_log <- v %in% c("Protein", "WBC", "RBC")
  
  if (use_log) {
    x <- log10(x_raw + 1)
  } else {
    x <- x_raw
  }
  
  df_plot <- data.frame(
    value = x,
    group = dat[[group_var]]
  ) %>%
    filter(!is.na(value), !is.na(group))
  

  n_per_group <- table(df_plot$group)
  if (any(n_per_group < 2)) {
    message("Skipping ", v, ": a group has fewer than two observations.")
    return(NULL)
  }
  

  # Final manuscript: all eight tests use the original-scale ranks.
  use <- complete.cases(dat[, c(v, group_var)])
  pval <- wilcox.test(dat[[v]][use] ~ dat[[group_var]][use],
                     exact = FALSE, correct = FALSE, alternative = "two.sided")$p.value
  p_stars <- if (pval < 0.001) {
    "***"
  } else if (pval < 0.01) {
    "**"
  } else if (pval < 0.05) {
    "*"
  } else {
    ""
  }
  

  p_tmp <- ggplot(
    df_plot,
    aes(x = group, y = value, fill = group)
  ) +
    geom_half_violin(trim = FALSE, side = "r")
  
  gb <- ggplot_build(p_tmp)
  y_violin_max <- max(gb$data[[1]]$y, na.rm = TRUE)
  y_violin_min <- min(gb$data[[1]]$y, na.rm = TRUE)
  y_range <- y_violin_max - y_violin_min
  

  line_y <- y_violin_max + 0.05 * y_range
  star_y <- y_violin_max + 0.10 * y_range
  

  ylim_lower <- y_violin_min - 0.2 * y_range
  ylim_upper <- y_violin_max + 0.2 * y_range
  

  xlim_lower <- 0.5
  xlim_upper <- 2.5
  line_x_start <- 1
  line_x_end   <- 2
  line_x_mid   <- 1.5
  

  plot_title <- if (use_log) paste0(v, " (log)") else v
  

  p <- ggplot(
    df_plot,
    aes(x = group, y = value, fill = group, color = group)
  ) +

    geom_half_violin(
      position = position_nudge(x = 0.1, y = 0),
      side = 'r',
      adjust = 1.2,
      trim = FALSE,
      color = NA,
      alpha = 0.6
    ) +

    geom_point(
      position = position_jitter(width = 0.15, height = 0),
      size = 1,
      alpha = 0.3
    ) +

    geom_boxplot(
      outlier.shape = NA,
      width = 0.1,
      alpha = 0.3
    ) +

    stat_boxplot(
      geom = "errorbar",
      width = 0.1,
      linewidth = 0.5
    ) +

    scale_color_manual(values = c("#1f77b4", "#d62728")) +
    scale_fill_manual(values = c("#1f77b4", "#d62728")) +

    {
      if (p_stars != "") {
        list(
          geom_segment(
            x     = line_x_start,
            xend  = line_x_end,
            y     = line_y,
            yend  = line_y,
            linewidth = 0.7,
            color     = "black"
          ),
          annotate(
            "text",
            x     = line_x_mid,
            y     = star_y,
            label = p_stars,
            size  = 5,
            vjust = 0

          )
        )
      } else {
        NULL
      }
    } +

    labs(
      title = plot_title,
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      panel.border = element_rect(linewidth = 0.8)
    ) +

    coord_cartesian(
      ylim   = c(ylim_lower, ylim_upper),
      xlim   = c(xlim_lower, xlim_upper),
      expand = FALSE
    )
  
  return(p)
}


plots<-lapply(csf_labels,function(v)make_violin_plot(v,data));names(plots)<-csf_labels
stopifnot(length(plots)==8L,!any(vapply(plots,is.null,logical(1))))
grobs<-lapply(plots,ggplot2::ggplotGrob)
# Equalize axis/title margins before the original grid.arrange-style assembly.
common_widths<-do.call(grid::unit.pmax,lapply(grobs,function(g)g$widths))
common_heights<-do.call(grid::unit.pmax,lapply(grobs,function(g)g$heights))
for(i in seq_along(grobs)) {
 grobs[[i]]$widths<-common_widths;grobs[[i]]$heights<-common_heights
 grobs[[i]]$layout$name[grobs[[i]]$layout$name=="panel"]<-paste0("data-panel-",i)
}
figure<-gridExtra::arrangeGrob(grobs=grobs,nrow=2,ncol=4)
# Measure the final R panel rectangles and block export if shared geometry differs.
source(file.path(project_dir,"tools/panel_alignment.R"),local=TRUE)
layout_path<-file.path(out,"figure2.alignment-layout.json")
write_gtable_panel_layout(figure,layout_path,width_in=11,height_in=8)
layout<-jsonlite::read_json(layout_path,simplifyVector=FALSE)
bbox<-do.call(rbind,lapply(layout$panels,function(p)unlist(p$bbox_pt)))
stopifnot(nrow(bbox)==8L,all(bbox[,3]>bbox[,1]),all(bbox[,4]>bbox[,2]))
stopifnot(diff(range(bbox[,3]-bbox[,1]))<=1.5,diff(range(bbox[,4]-bbox[,2]))<=1.5)
prefix<-file.path(out,"figure2")
draw_figure<-function(){grid::grid.newpage();grid::grid.draw(figure)}
pdf(paste0(prefix,".pdf"),width=11,height=8,family="Helvetica",useDingbats=FALSE);draw_figure();dev.off()
png(paste0(prefix,".png"),width=11,height=8,units="in",res=300,type="cairo",bg="white");draw_figure();dev.off()
tiff(paste0(prefix,".tiff"),width=11,height=8,units="in",res=600,type="cairo",compression="lzw",bg="white");draw_figure();dev.off()
svglite::svglite(paste0(prefix,".svg"),width=11,height=8,bg="white");draw_figure();dev.off()
stats<-comparison_table(d,csf_vars,csf_labels,FALSE,TRUE);stats$Display_transform<-ifelse(csf_labels%in%c("Protein","WBC","RBC"),"log10(x+1)","none")
stats$Stars<-ifelse(stats$P_value<.001,"***",ifelse(stats$P_value<.01,"**",ifelse(stats$P_value<.05,"*","")))
stopifnot(stats$Stars[stats$Variable=="Chloride"]=="")
write_csv(stats,file.path(out,"figure2_source_statistics.csv"));write_csv(data,file.path(out,"figure2_source_data.csv"))
writeLines(c("Eight panels retain the provided half-violin, scatter, boxplot and colors.","All observations with available values are shown. No sampling. Jitter seed: 1963580.","The final manuscript uses Wilcoxon tests for all eight comparisons; the earlier template's t-test branch is not used.","Only protein, WBC and RBC are log10(x+1)-transformed for display. Units refer to the original measurements."),file.path(out,"figure2_notes.txt"))
