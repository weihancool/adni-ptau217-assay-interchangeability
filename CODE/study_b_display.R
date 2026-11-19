sb_display_palette <- c(ink = "#253444", blue = "#254F70", red = "#B44145",
                        adjacent = "#8597A6", quiet = "#BBC5CC", grid = "#E7EBEE")
sb_theme_classic <- function(base_size = 8.5, base_family = "sans", ...) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family, ...) +
    ggplot2::theme(text = ggplot2::element_text(colour = sb_display_palette[["ink"]]),
                   axis.text = ggplot2::element_text(size = max(7,base_size - 1),colour = sb_display_palette[["ink"]]),
                   axis.title = ggplot2::element_text(size = base_size),
                   axis.line = ggplot2::element_line(colour = "#74818D",linewidth = .3),
                   axis.ticks = ggplot2::element_line(colour = "#74818D",linewidth = .3),
                   axis.ticks.length = grid::unit(1.1,"mm"),
                   plot.title = ggplot2::element_text(size = 11,face = "bold",margin = ggplot2::margin(b = 4)),
                   plot.title.position = "plot",plot.margin = ggplot2::margin(4,6,5,5),
                   panel.grid = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white",colour = NA))
}
sb_theme_minimal <- function(base_size = 8.5, base_family = "sans", ...) {
  sb_theme_classic(base_size,base_family,...) +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = sb_display_palette[["grid"]],linewidth = .3))
}
sb_ggsave <- function(filename,plot = ggplot2::last_plot(),device = NULL,
                      path = NULL,scale = 1,width = NA,height = NA,units = "in",
                      dpi = 600,limitsize = TRUE,bg = "white",...) {
  if (exists("sb_assert_code_current",mode = "function")) sb_assert_code_current()
  output_path <- if (is.null(path)) filename else file.path(path,filename)
  if (exists(".sb_reporting_state",inherits=TRUE) &&
      normalizePath(output_path,winslash="/",mustWork=FALSE) %in% names(.sb_reporting_state$inputs))
    stop("Refusing to overwrite a consumed frozen figure: ",output_path,call.=FALSE)
  ext <- tolower(tools::file_ext(filename))
  if (is.null(device)) device <- switch(ext,
                                        pdf = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf,
                                        png = grDevices::png,tif = grDevices::tiff,tiff = grDevices::tiff,NULL)
  if (ext == "png" && capabilities("cairo") && identical(device,grDevices::png))
    device <- function(filename,width,height,bg,...) {
      grDevices::png(filename = filename,width = width,height = height,units = "in",
                     res = dpi,bg = bg,type = "cairo-png",...)
    }
  if (is.numeric(dpi)) dpi <- max(600,dpi)
  ggplot2::ggsave(filename,plot = plot,device = device,path = path,scale = scale,
                  width = width,height = height,units = units,dpi = dpi,limitsize = limitsize,bg = bg,...)
  full_path <- if (is.null(path)) filename else file.path(path,filename)
  if (!file.exists(full_path) || is.na(file.info(full_path)$size) || file.info(full_path)$size <= 0)
    stop("Figure export absent/empty: ",full_path,call. = FALSE)
  if (ext == "png") {
    con <- file(full_path,"rb")
    dims <- tryCatch({readBin(con,"raw",n = 16);readBin(con,"integer",n = 2,size = 4,endian = "big")},
                     finally = close(con))
    if (length(dims) != 2L || any(dims < 100L))
      stop("PNG dimensions are implausibly small; check device units: ",full_path,call. = FALSE)
  }
  if (exists("sb_stamp_result",mode = "function")) sb_stamp_result(full_path)
  invisible(full_path)
}
sb_export_figure <- function(plot,base,width = 190,height = 155) {
  for (ext in c("pdf","png","tiff")) {
    args <- list(filename = paste0(base,".",ext),plot = plot,width = width,height = height,units = "mm")
    if (ext == "tiff") args$compression <- "lzw"
    do.call(sb_ggsave,args)
  }
  invisible(base)
}
sb_require_results <- function(files,producer,res = study_b_dir) {
  paths <- file.path(res,files)
  csv_ok <- file.exists(paths); provenance_ok <- file.exists(paste0(paths,".provenance.rds"))
  if (any(!csv_ok | !provenance_ok)) {
    missing <- paste0(files[!csv_ok | !provenance_ok]," [CSV=",csv_ok[!csv_ok | !provenance_ok],
                      ", provenance=",provenance_ok[!csv_ok | !provenance_ok],"]")
    stop("Required results are incomplete in:\n",normalizePath(res,winslash = "/",mustWork = FALSE),
         "\n",paste(missing,collapse = "\n"),"\nProducer: ",producer,
         "\nIf this is the wrong project, correct STUDY_B_ROOT and run 01_config.R again.",
         "\nIf it is the correct project, run the producer after its prerequisites.",
         "\nA CSV alone does not satisfy the frozen-result contract.",call. = FALSE)
  }
  invisible(paths)
}
sb_display_map <- function() {
  data.frame(
    figure = c(rep("Figure 1",4),rep("Figure 2",3),rep("Figure 3",4)),
    panel = c("A","B","C","D","A","B","C","A","B","C","D"),
    evidence_key = c("PET-discrimination-by-assay","joint-min-max-Youden-margin",
                     "joint-four-assay-binary-disagreement","joint-three-level-composition",
                     "same-draw-pairwise-binary-switching","all-occasion-margin-stratified-switching-with-residual-benchmark",
                     "same-draw-switching-by-population","longitudinal-change-joint-distribution",
                     "longitudinal-change-magnitude-ecdf","longitudinal-pair-reversal-and-dominance",
                     "longitudinal-category-controls-and-paired-excess"),
    population = c("PET-matched primary anchors",rep("Complete primary anchors",3),
                   "Complete primary anchors","All eligible collection occasions",
                   "Primary anchors: full cohort and MCI/dementia",
                   rep("Aligned first/last draws in repeated participants",4)),
    stringsAsFactors = FALSE)
}
sb_finalize_display_outputs <- function(res = study_b_dir) {
  all_files <- list.files(res,pattern = "^Table_S[0-9]+_.*[.]csv$",full.names = FALSE)
  all_files <- all_files[!grepl("_final[.]csv$|_display[.]csv$|bootstrap_replicates",all_files)]
  canonical <- c("Table_S1_assay_decision_prevalence.csv","Table_S2_pairwise_agreement_full.csv",
                 "Table_S3_gray_zone_patterns_full.csv","Table_S4_variance_components_full.csv",
                 "Table_S5_fixed_thresholds.csv","Table_S6_sensitivity_window_cutoff_full.csv",
                 "Table_S7_A1_assay_switching_consequences.csv","Table_S7_A1_residual_variance_benchmark.csv",
                 "Table_S8_A2_symptomatic_switching_consequences.csv",
                 "Table_S9_A3_longitudinal_pseudo_switch.csv","Table_S10_B3_hypothetical_prevalence_by_assay.csv",
                 "Table_S11_B2_variance_explanation.csv")
  sb_require_results(canonical,"current supplementary-table producers in scripts 05-30 (see README)",res)
  source_data <- list(); excluded <- list()
  for (f in all_files) {
    z <- tryCatch(sb_read_result(file.path(res,f),check.names = FALSE),error = function(e) {
      if (f %in% canonical) stop(conditionMessage(e),call. = FALSE)
      excluded[[length(excluded)+1L]] <<- data.frame(source_file = f,reason = conditionMessage(e))
      NULL
    })
    if (!is.null(z)) source_data[[f]] <- z
  }
  excluded_df <- if (length(excluded)) dplyr::bind_rows(excluded) else
    data.frame(source_file = character(),reason = character())
  sb_save_csv(excluded_df,file.path(res,"EXCLUDED_STALE_SUPPLEMENT_SOURCES.csv"),row.names = FALSE,na = "")
  if (nrow(excluded_df)) warning(nrow(excluded_df)," old/optional supplementary source(s) excluded; see EXCLUDED_STALE_SUPPLEMENT_SOURCES.csv.",call. = FALSE)
  all_files <- names(source_data)
  index <- list()
  for (i in 1:11) {
    files <- all_files[grepl(paste0("^Table_S",i,"_"),all_files)]
    if (!length(files)) stop("No source files for supplementary Table S",i,". Complete scripts 01-30.",call. = FALSE)
    parts <- lapply(files,function(f) {
      z <- source_data[[f]]
      z$source_file <- f; z$table_section <- sub(paste0("^Table_S",i,"_"),"",sub("[.]csv$","",f))
      z$section_role <- if (grepl("audit|diagnostic|contract|coverage|missing",f,ignore.case = TRUE)) "Audit/support" else "Analysis"
      z[] <- lapply(z,as.character)
      index[[length(index)+1L]] <<- data.frame(table = paste0("Table S",i),source_file = f,
                                               rows = nrow(z),section_role = if (nrow(z)) z$section_role[1] else "Empty section",
                                               stringsAsFactors = FALSE)
      z
    })
    combined <- dplyr::bind_rows(parts)
    sb_save_csv(combined,file.path(res,paste0("Table_S",i,"_final.csv")),row.names = FALSE,na = "")
  }
  sb_save_csv(dplyr::bind_rows(index),file.path(res,"SUPPLEMENT_TABLE_INDEX.csv"),row.names = FALSE,na = "")
  map <- sb_display_map()
  stopifnot(!anyDuplicated(map[c("figure","panel")]),!anyDuplicated(map$evidence_key))
  sb_save_csv(map,file.path(res,"MANUSCRIPT_DISPLAY_MAP.csv"),row.names = FALSE,na = "")
  stems <- c("Figure_1_interchangeability","Figure_2_switching_consequences",
             "Figure_3_longitudinal_switch_consequences","Figure_S1_flowchart",
             "Figure_S2_gray_zone_patterns_top15","Figure_S3_threshold_distributions",
             "Figure_S4_sensitivity_ratio","Figure_S5_variance_decomposition",
             "Figure_S6_longitudinal_persistence","Figure_S7_remeasurement_vs_switching",
             "Figure_S8_A3_small_change_sensitivity")
  valid_export <- function(file) {
    path <- file.path(res,file)
    tryCatch({
      meta <- sb_validate_result(path)
      if (startsWith(file,"Figure_2_switching_consequences.") &&
          !identical(meta$producer_script,"16_figure2_switching_consequences.R"))
        return(FALSE)
      .sb_note_input(path)
      file.info(path)$size > 0
    }, error = function(e) FALSE)
  }
  present <- vapply(stems,function(stem) {
    files <- paste0(stem,c(".pdf",".png",".tiff"))
    all(vapply(files,valid_export,logical(1)))
  },logical(1))
  audit <- data.frame(display = c("Figure 1","Figure 2","Figure 3",paste0("Figure S",1:8)),
                      all_three_formats_present = present)
  sb_save_csv(audit,file.path(res,"FINAL_DISPLAY_EXPORT_CHECKS.csv"),row.names = FALSE,na = "")
  if (any(!present)) stop("Missing, stale or unverified figure exports: ",paste(audit$display[!present],collapse = ", "),call. = FALSE)
  invisible(audit)
}
