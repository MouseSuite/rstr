# Rodent Statistics Toolbox in R (rstr)
# Copyright (C) 2025 The Regents of the University of California
# Creator: Shantanu H. Joshi, Department of Neurology, Ahmanson Lovelace Brain Mapping Center, UCLA
#
# This program is free software; you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation; version 2.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License version 2 for more details.
#
# You should have received a copy of the GNU General Public License along with this program;
# if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

#' Check that subject directory and demographics files exist
#' @param object object of type \code{RstrOutput}
#'
check_files <- function(object){
  if (!dir.exists(object@subjdir)) {
    stop(sprintf("Subjects directory %s does not exist.\n", object@subjdir), call. = FALSE)
  }

  if (!file.exists(object@csv)) {
    stop(sprintf("Demographics csv file %s does not exist.\n", object@csv), call. = FALSE)
  }
}

#' S4 class for saving results of statistical analysis
#' @slot outdir output directory to save the results
#'
RstrOutput <- setClass(
  "RstrOutput",
  slots = list(
    outdir = "character"
  ),
  validity = check_files
)

setMethod("initialize", valueClass = "RstrOutput", signature = "RstrOutput", function(.Object, outdir) {
  if (!dir.exists(outdir)) {
    dir.create(outdir)
    .Object@outdir <- outdir
  }
  else {
    .Object@outdir <- outdir
    message(sprintf("The output directory %s already exists.", outdir))
  }
  return(.Object)
})

#' Generic save function for \code{RstrOutput}
#' @param rstr_out object of type \code{RstrOutput}
#' @param rstr_data object of type \code{RstrData}
#' @param rstr_model object of type \code{RstrModel}
#' @param overwrite logical parameter denoting if existing output directory should be overwritten or not (default is FALSE)
#' @param ... Extra named arguments passed to save_out
#' @details
#' For the most part, the user will never have to call this function directly.
#' Instead the user should call \code{\link{save_rstr_out}}.
#' @seealso \code{\link{save_rstr_out}}
#' @export
setGeneric("save_out", valueClass = "RstrOutput", function(rstr_out, rstr_data, rstr_model, overwrite = FALSE, ...) {
  standardGeneric("save_out")
})

RstrSBAOutput <- setClass(
  "RstrSBAOutput",
  contains = "RstrOutput"
)

RstrTBMOutput <- setClass(
  "RstrTBMOutput",
  contains = "RstrOutput"
)

RstrDBAOutput <- setClass(
  "RstrDBAOutput",
  contains = "RstrOutput"
)

RstrROIOutput <- setClass(
  "RstrROIOutput",
  contains = "RstrOutput"
)


#' @rdname save_out
#' @inheritParams save_out
setMethod("save_out", valueClass = "RstrSBAOutput", signature = "RstrSBAOutput", function(rstr_out, rstr_data, rstr_model, overwrite = F) {

  # If output directory is not empty, then empty if overwrite is true or stop if overwrite is false
  if (overwrite == TRUE){
    delete_and_recreate_dir(rstr_out@outdir)
  }
  # else {
  #   if (length(list.files(rstr_out@outdir, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) != 0){
  #     stop(sprintf("Output directory %s is not empty.\n", rstr_out@outdir), call. = FALSE)
  #   }
  # }

  log_pvalues <- log10_transform(rstr_model@pvalues)
  outdir <- rstr_out@outdir
  log_pvalues_adjusted <- log10_transform(rstr_model@pvalues_adjusted)
  rstr_model@tvalues[abs(log_pvalues) <= -1*log10(0.05)] <- 0

  switch(rstr_model@model_type,
         rstr_anova = {
           rstr_cmap <- save_rstr_color_files(log_pvalues, rstr_model@main_effect, "log_pvalues", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(log_pvalues, rstr_model@main_effect, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(log_pvalues_adjusted, rstr_model@main_effect, "log_pvalues_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(log_pvalues_adjusted, rstr_model@main_effect, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@tvalues, rstr_model@main_effect, "tvalues", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@tvalues, rstr_model@main_effect, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@tvalues_adjusted, rstr_model@main_effect, "tvalues_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@tvalues_adjusted, rstr_model@main_effect, rstr_cmap, rstr_data, rstr_model, outdir)
           save_rstr_rds(rstr_model@pvalues, rstr_model@main_effect, "pvalues", rstr_data, rstr_model, outdir) # Save pvalues as a rds file
         },
         rstr_lm = {
           rstr_cmap <- save_rstr_color_files(log_pvalues, rstr_model@main_effect, "log_pvalues", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(log_pvalues, rstr_model@main_effect, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(log_pvalues_adjusted, rstr_model@main_effect, "log_pvalues_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(log_pvalues_adjusted, rstr_model@main_effect, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@tvalues, rstr_model@main_effect, "tvalues", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@tvalues, rstr_model@main_effect, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@tvalues_adjusted, rstr_model@main_effect, "tvalues_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@tvalues_adjusted, rstr_model@main_effect, rstr_cmap, rstr_data, rstr_model, outdir)
           save_rstr_rds(rstr_model@pvalues, rstr_model@main_effect, "pvalues", rstr_data, rstr_model, outdir) # Save pvalues as a rds file
           },
         rstr_corr = {
           rstr_model@corr_values[abs(log_pvalues) <= -1*log10(0.05)] <- 0
           rstr_cmap <- save_rstr_color_files(log_pvalues, rstr_model@corr_var, "log_pvalues", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(log_pvalues, rstr_model@corr_var, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(log_pvalues_adjusted, rstr_model@corr_var, "log_pvalues_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(log_pvalues_adjusted, rstr_model@corr_var, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@tvalues, rstr_model@corr_var, "tvalues", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@tvalues, rstr_model@corr_var, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@tvalues_adjusted, rstr_model@corr_var, "tvalues_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@tvalues_adjusted, rstr_model@corr_var, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@corr_values, rstr_model@corr_var, "corr_values", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@corr_values, rstr_model@corr_var, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@corr_values_masked_adjusted, rstr_model@corr_var, "corr_values_masked_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@corr_values_masked_adjusted, rstr_model@corr_var, rstr_cmap, rstr_data, rstr_model, outdir)
           save_rstr_rds(rstr_model@pvalues, rstr_model@corr_var, "pvalues", rstr_data, rstr_model, outdir) # Save pvalues as a rds file
           },
         pairedttest = {
           rstr_cmap <- save_rstr_color_files(log_pvalues, rstr_model@group_var, "log_pvalues", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(log_pvalues, rstr_model@group_var, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(log_pvalues_adjusted, rstr_model@group_var, "log_pvalues_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(log_pvalues_adjusted, rstr_model@group_var, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@tvalues, rstr_model@group_var, "tvalues", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@tvalues, rstr_model@group_var, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@tvalues_adjusted, rstr_model@group_var, "tvalues_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@tvalues_adjusted, rstr_model@group_var, rstr_cmap, rstr_data, rstr_model, outdir)
           save_rstr_rds(rstr_model@pvalues, rstr_model@group_var, "pvalues", rstr_data, rstr_model, outdir) # Save pvalues as a rds file
         },
         unpairedttest = {
           rstr_cmap <- save_rstr_color_files(log_pvalues, rstr_model@group_var, "log_pvalues", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(log_pvalues, rstr_model@group_var, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(log_pvalues_adjusted, rstr_model@group_var, "log_pvalues_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(log_pvalues_adjusted, rstr_model@group_var, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@tvalues, rstr_model@group_var, "tvalues", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@tvalues, rstr_model@group_var, rstr_cmap, rstr_data, rstr_model, outdir)
           rstr_cmap <- save_rstr_color_files(rstr_model@tvalues_adjusted, rstr_model@group_var, "tvalues_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_surface_both_hemi(rstr_model@tvalues_adjusted, rstr_model@group_var, rstr_cmap, rstr_data, rstr_model, outdir)
           save_rstr_rds(rstr_model@pvalues, rstr_model@group_var, "pvalues", rstr_data, rstr_model, outdir) # Save pvalues as a rds file
         }
  )

  # Copy modelspec file to the output directory
  file.copy(rstr_model@mspec_file, rstr_out@outdir)
  return(rstr_out)
  }
)


#' @rdname save_out
#' @inheritParams save_out
#' @param nclusters numeric parameter denoting number of clusters (default is 10)
setMethod("save_out", valueClass = "RstrTBMOutput", signature = "RstrTBMOutput", function(rstr_out, rstr_data, rstr_model, overwrite = F, nclusters = 10) {

  # If output directory is not empty, then empty if overwrite is true or stop if overwrite is false
  if (overwrite == TRUE){
    delete_and_recreate_dir(rstr_out@outdir)
  } else {
    if (length(list.files(rstr_out@outdir, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) != 0){
      stop(sprintf("Output directory %s is not empty.\n", rstr_out@outdir), call. = FALSE)
    }
  }

  log_pvalues <- rep(1, length(rstr_data@atlas_image))
  log_pvalues[rstr_data@mask_idx] <- log10_transform(rstr_model@pvalues)
  dim(log_pvalues) <- dim(rstr_data@atlas_image)

  log_pvalues_adjusted <- rep(1, length(rstr_data@atlas_image))
  log_pvalues_adjusted[rstr_data@mask_idx] <- log10_transform(rstr_model@pvalues_adjusted)
  dim(log_pvalues_adjusted) <- dim(rstr_data@atlas_image)
  outdir <- rstr_out@outdir

  tvalues <- rep(0, length(rstr_data@atlas_image))
  tvalues[rstr_data@mask_idx] <- rstr_model@tvalues
  dim(tvalues) <- dim(rstr_data@atlas_image)

  tvalues_adjusted <- rep(0, length(rstr_data@atlas_image))
  tvalues_adjusted[rstr_data@mask_idx] <- rstr_model@tvalues_adjusted
  dim(tvalues_adjusted) <- dim(rstr_data@atlas_image)
  measure <- NULL
  switch(rstr_model@model_type,
         rstr_anova = {
           var_name = rstr_model@main_effect
         },
         rstr_lm = {
           var_name = rstr_model@main_effect
         },
         rstr_corr = {
           corr_values <- rep(0, length(rstr_data@atlas_image))
           corr_values[rstr_data@mask_idx] <- rstr_model@corr_values
           dim(corr_values) <- dim(rstr_data@atlas_image)

           corr_values_masked_adjusted <- rep(0, length(rstr_data@atlas_image))
           corr_values_masked_adjusted[rstr_data@mask_idx] <- rstr_model@corr_values_masked_adjusted
           dim(corr_values_masked_adjusted) <- dim(rstr_data@atlas_image)

           var_name = rstr_model@corr_var
         },
         pairedttest = {
           var_name = rstr_model@group_var
         },
         unpairedttest = {
           var_name = rstr_model@group_var
         }
  )

  if (rstr_model@model_type == "rstr_corr") {
    save_vol_stats_out(log_pvalues, log_pvalues_adjusted, tvalues, tvalues_adjusted,
                       var_name, rstr_data, rstr_model, outdir, corr_values, corr_values_masked_adjusted)
  } else {
    save_vol_stats_out(log_pvalues, log_pvalues_adjusted, tvalues, tvalues_adjusted, var_name, rstr_data, rstr_model, outdir)
  }

  voxelcoord <- get_voxelcoord(rstr_out, rstr_data, rstr_model, outdir, nclusters)

  # Check if voxelcoord is empty
  if (length(voxelcoord) == 0 | voxelcoord[[1]][1] == -1) {
    sink(file.path(outdir,  sprintf("report_%s_%s.Rmd", rstr_model@model_type, var_name)), type = "output")
    cat("---\n")
    cat("title: RSTR Report\n")
    cat("output: html_document\n")
    cat("---\n\n\n")
    cat("No detected clusters above significance threshold.")
    sink()
    rmarkdown::render(file.path(outdir,  sprintf("report_%s_%s.Rmd", rstr_model@model_type, var_name)))
    stop("No detected clusters above significance threshold.")
  }

  # add create an R6 class function from here
  rstrmd_volout <- RstrRmdVolumeOutput$new()
  rstrmd_volout$save_out(rstr_data, rstr_model, voxelcoord = voxelcoord, outdir)


    # Copy modelspec file to the output directory
  file.copy(rstr_model@mspec_file, rstr_out@outdir)
  invisible(rstr_out)
  }
)


#' @rdname save_out
#' @inheritParams save_out
setMethod("save_out", valueClass = "RstrDBAOutput", signature = "RstrDBAOutput", function(rstr_out, rstr_data, rstr_model, overwrite = F, nclusters = 10) {

  # If output directory is not empty, then empty if overwrite is true or stop if overwrite is false
  if (overwrite == TRUE){
    delete_and_recreate_dir(rstr_out@outdir)
  } else {
    if (length(list.files(rstr_out@outdir, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) != 0){
      stop(sprintf("Output directory %s is not empty.\n", rstr_out@outdir), call. = FALSE)
    }
  }

  log_pvalues <- rep(1, length(rstr_data@atlas_image))
  log_pvalues[rstr_data@mask_idx] <- log10_transform(rstr_model@pvalues)
  dim(log_pvalues) <- dim(rstr_data@atlas_image)

  log_pvalues_adjusted <- rep(1, length(rstr_data@atlas_image))
  log_pvalues_adjusted[rstr_data@mask_idx] <- log10_transform(rstr_model@pvalues_adjusted)
  dim(log_pvalues_adjusted) <- dim(rstr_data@atlas_image)
  outdir <- rstr_out@outdir

  tvalues <- rep(0, length(rstr_data@atlas_image))
  tvalues[rstr_data@mask_idx] <- rstr_model@tvalues
  dim(tvalues) <- dim(rstr_data@atlas_image)

  tvalues_adjusted <- rep(0, length(rstr_data@atlas_image))
  tvalues_adjusted[rstr_data@mask_idx] <- rstr_model@tvalues_adjusted
  dim(tvalues_adjusted) <- dim(rstr_data@atlas_image)
  measure <- NULL
  switch(rstr_model@model_type,
         rstr_anova = {
           var_name = rstr_model@main_effect
         },
         rstr_lm = {
           var_name = rstr_model@main_effect
         },
         rstr_corr = {
           corr_values <- rep(0, length(rstr_data@atlas_image))
           corr_values[rstr_data@mask_idx] <- rstr_model@corr_values
           dim(corr_values) <- dim(rstr_data@atlas_image)

           corr_values_masked_adjusted <- rep(0, length(rstr_data@atlas_image))
           corr_values_masked_adjusted[rstr_data@mask_idx] <- rstr_model@corr_values_masked_adjusted
           dim(corr_values_masked_adjusted) <- dim(rstr_data@atlas_image)

           var_name = rstr_model@corr_var
         },
         pairedttest = {
           var_name = rstr_model@group_var
         },
         unpairedttest = {
           var_name = rstr_model@group_var
         }
  )

  if (rstr_model@model_type == "rstr_corr") {
    save_vol_stats_out(log_pvalues, log_pvalues_adjusted, tvalues, tvalues_adjusted,
                       var_name, rstr_data, rstr_model, outdir, corr_values, corr_values_masked_adjusted)
  } else {
    save_vol_stats_out(log_pvalues, log_pvalues_adjusted, tvalues, tvalues_adjusted, var_name, rstr_data, rstr_model, outdir)
  }

  voxelcoord <- get_voxelcoord(rstr_out, rstr_data, rstr_model, outdir, nclusters)

  # Check if voxelcoord is empty
  if (length(voxelcoord) == 0 | voxelcoord[[1]][1] == -1) {
    sink(file.path(outdir,  sprintf("report_%s_%s.Rmd", rstr_model@model_type, var_name)), type = "output")
    cat("---\n")
    cat("title: RSTR Report\n")
    cat("output: html_document\n")
    cat("---\n\n\n")
    cat("No detected clusters above significance threshold.")
    sink()
    rmarkdown::render(file.path(outdir,  sprintf("report_%s_%s.Rmd", rstr_model@model_type, var_name)))
    stop("No detected clusters above significance threshold.")
  }

  # Create a new R6 class object here
  rstrmd_volout <- RstrRmdVolumeOutput$new()
  rstrmd_volout$save_out(rstr_data, rstr_model, voxelcoord = voxelcoord, outdir)

  # Copy modelspec file to the output directory
  file.copy(rstr_model@mspec_file, rstr_out@outdir)
  invisible(rstr_out)
}
)

#' @rdname save_out
#' @inheritParams save_out
setMethod("save_out", valueClass = "RstrROIOutput", signature = "RstrROIOutput", function(rstr_out, rstr_data, rstr_model, overwrite = F) {

  # # Create the output directory
  # if (is.null(outdir)) {
  #   # Outputted directory (if not specified) is of the same format (csv or tsv) as the inputted file
  #   if (tools::file_ext(file.path(csv)) == "csv"){
  #     outdir <- paste0(tools::file_path_sans_ext(csv), "_roidata.csv")
  #     data_separator = ','
  #   } else if (tools::file_ext(file.path(csv)) == "tsv"){
  #     outdir <- paste0(tools::file_path_sans_ext(csv), "_roidata.tsv")
  #     data_separator = '\t'
  #   }
  #   cat(sprintf('Output directory is not specified. Using %s to save outputs.\n', outdir))
  # }
  # else {
  #   dir.create(file.path(outdir), showWarnings = FALSE)
  # }
  #
  # write.table(rstr_data@demographics, outdir, sep = data_separator,row.names = FALSE)

  # If output directory is not empty, then empty if overwrite is true or stop if overwrite is false
  if (overwrite == TRUE){
    delete_and_recreate_dir(rstr_out@outdir)
  } else {
    if (length(list.files(rstr_out@outdir, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) != 0){
      stop(sprintf("Output directory %s is not empty.\n", rstr_out@outdir), call. = FALSE)
    }
  }

  # Get the absolute path of outdir
  outdir <- tools::file_path_as_absolute(rstr_out@outdir)

  # Copy demographics csv file to output directory
  csvfilename <- paste0(file.path(rstr_out@outdir, tools::file_path_sans_ext(basename(rstr_data@csv))),".csv")
  write.csv(rstr_data@demographics, csvfilename)

  nb_header <- "---\ntitle: 'MouseSuite ROI statistical analysis report'\noutput: html_document\n---"

  nb_libraries <-"```{r librar_cmds, echo=FALSE}\n"
  nb_libraries <- paste(nb_libraries, "\nlibrary('rstr')\nlibrary('ggplot2')\n", sep="")
  nb_libraries <- paste(nb_libraries, "\n```\n", sep="")

  nb_data_header_one <- "The following command loads the data."

  # Add commands to load the data
  nb_data_command_one <- sprintf("\n```{r message=FALSE, warning=FALSE, results='hide', data_command_one}\n")
  nb_data_command_one <- paste(nb_data_command_one,"\n",
                            rstr_data@load_data_command,"\n",sep = "")
  nb_data_command_one <- paste(nb_data_command_one, "\n```\n", sep = "")


  nb_load_data <- sprintf("```{r echo=FALSE, message=FALSE, warning=FALSE, load_data}\n")
  nb_load_data <- paste0(nb_load_data, "\nDT::datatable(rstr_data@demographics[,-(ncol(rstr_data@demographics))], rownames = FALSE)\n")
  nb_load_data <- paste0(nb_load_data, "\n```\n")


  nb_data_header_two <- "The following command creates the model used to analyze the data."

  nb_data_command_two <- sprintf("\n```{r message=FALSE, warning=FALSE, results='hide', data_command_two}\n")
  nb_data_command_two <- paste(nb_data_command_two, "\n",rstr_model@load_data_command,"\n")
  nb_data_command_two <- paste0(nb_data_command_two, "\n```\n")

  nb_data_header_three <- "The final command (below) was used to render this document."

  command_save_rstr_out <- sprintf("save_rstr_out(rstr_data, rstr_model, outdir = '%s')", rstr_out@outdir)

  nb_data_command_three <-  sprintf("\n```{r eval=FALSE, message=FALSE, warning=FALSE, data_command_three}\n")
  nb_data_command_three <- paste(nb_data_command_three, "\n", command_save_rstr_out,"\n")
  nb_data_command_three <- paste0(nb_data_command_three, "\n```\n")


  nb_commands <- lapply(as.character(rstr_data@roiids),function(x){sprintf("```{r warning=FALSE, run_command_%s}\n",x)})
  nb_plots <- lapply(as.character(rstr_data@roiids),function(x){sprintf("```{r echo=FALSE, message=FALSE, warning=FALSE, plot_%s}\n",x)})
  if (rstr_model@model_type != 'rstr_corr'){
    nb_calculations <- lapply(as.character(rstr_data@roiids),function(x){sprintf("```{r echo=FALSE, message=FALSE, warning=FALSE, pval_%s}\n",x)})
  }

  selected_col <- rep(NA, length(rstr_data@roiids))
  for (i in 1:length(rstr_data@roiids)){
    selected_col[i] <- paste0(as.character(get_roi_tag(read_label_desc(rstr_data@roilabeldescfile),rstr_data@roiids[i])), "(",rstr_data@roiids[i],")")
    rstr_data@demographics[,selected_col[i]]
  }

  if (rstr_model@model_type != 'pairedttest' & rstr_model@model_type != 'unpairedttest') {
    if (rstr_model@model_type=="rstr_lmer"){ comparison_stat = "Chisq"}
    else {comparison_stat = "F"}

    for (m in 1:length(rstr_data@roiids)){
      if (rstr_model@model_type != 'rstr_corr'){
        for (i in 1:length(rstr_model@stats_commands[[m]])) {
          nb_commands[[m]] <- paste0(nb_commands[[m]], rstr_model@stats_commands[[m]][i], "\n")
        }
        nb_commands[[m]] <- paste0(nb_commands[[m]], "```\n\n")
        nb_calculations[[m]] <- paste0(nb_calculations[[m]],"anova_table <- ",
                                       substr(rstr_model@stats_commands[[m]][3],nchar("pander::pander(")+1,nchar(rstr_model@stats_commands[[m]][3])-1),
                                       "\np_val_",rstr_data@roiids[m],"<- round(anova_table$`Pr(>",comparison_stat,")`[2],digits=4)\n",
                                       "pval_string <- paste('pvalue:', as.character(p_val_",rstr_data@roiids[m],"))\n```\n\n")
      }
      if (class(rstr_data@demographics[,gsub("([[:alnum:]_]+).*", "\\1", rstr_model@fullmodel)])=="integer"|class(rstr_data@demographics[,gsub("([[:alnum:]_]+).*", "\\1", rstr_model@fullmodel)])=="double"|rstr_model@model_type == 'rstr_corr'){
        if (rstr_model@model_type == "rstr_corr"){
          x_var = rstr_model@corr_var
          annotate_label = paste0("paste('corr val:',",round(rstr_model@corr_values[m],5),")")
        } else {
          x_var = gsub('([[:alnum:]_]+).*', '\\1', rstr_model@fullmodel)
          annotate_label = paste0("paste('pvalue:', as.character(p_val_",rstr_data@roiids[m],"))")
        }
        nb_plots[[m]]<-paste0(nb_plots[[m]],"ggplot2::ggplot(data=rstr_data@demographics, ggplot2::aes(x=",
                              x_var,", y = `",selected_col[m],"`)) + ggplot2::geom_point() + ggplot2::ggtitle('",
                              as.character(get_roi_name(label_desc_df = read_label_desc(rstr_data@roilabeldescfile),roiid=rstr_data@roiids[m])[[1]]),
                              " ", rstr_data@roimeas,
                              " vs ", x_var, "') + ggplot2::labs(y='",
                              as.character(get_roi_name(label_desc_df = read_label_desc(rstr_data@roilabeldescfile),roiid=rstr_data@roiids[m])[[1]]),
                              "') + ggplot2::geom_smooth(method=lm, se=TRUE) +
                              ggplot2::xlim(",min(rstr_data@demographics[,x_var]),", ",max(rstr_data@demographics[,x_var])+6,") +
                              ggplot2::theme(axis.title=ggplot2::element_text(size=16,face='bold')) +
                              ggplot2::theme(plot.title=ggplot2::element_text(size=18,face='bold')) +
                              ggplot2::annotate('label',x=",max(rstr_data@demographics[,x_var])+2,",y= max(rstr_data@demographics$`",selected_col[m],"`)",
                              ",label= ",annotate_label,")\n
                              ggplot2::ggsave(filename='",
                              paste0(as.character(get_roi_tag(read_label_desc(rstr_data@roilabeldescfile),rstr_data@roiids[m]))), "_roi",rstr_data@roiids[m],
                              "_", rstr_data@roimeas,
                              "_vs_", x_var, ".pdf',device='pdf')\n```\n")
      } else if (class(rstr_data@demographics[,gsub("([[:alnum:]_]+).*", "\\1", rstr_model@fullmodel)])=="factor"){
        nb_plots[[m]]<-paste0(nb_plots[[m]],
                              "mean_lengths <- rep(NA, length(levels(rstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', rstr_model@fullmodel),")))
                              std_devs <- rep(NA, length(levels(rstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', rstr_model@fullmodel),")))
                              for (i in 1:length(levels(rstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', rstr_model@fullmodel),"))){
                              mean_lengths[i] <-  mean(rstr_data@demographics[rstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', rstr_model@fullmodel),
                              " == levels(rstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', rstr_model@fullmodel),")[i],]$`",
                              selected_col[m],"`)
                              std_devs[i] <- sd (rstr_data@demographics[rstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', rstr_model@fullmodel),
                              " == levels(rstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', rstr_model@fullmodel),")[i],]$`",
                              selected_col[m],"`)
                              } \n",
"modified_df <- data.frame(len = mean_lengths, sd = std_devs,",gsub('([[:alnum:]_]+).*', '\\1', rstr_model@fullmodel)," = levels(rstr_data@demographics$",
gsub('([A-Z_a-z]+).*', '\\1', rstr_model@fullmodel),"))\n",
"ggplot2::ggplot(data=modified_df, ggplot2::aes(x = ",gsub('([[:alnum:]_]+).*', '\\1', rstr_model@fullmodel),
", y = len, color = ",gsub('([[:alnum:]_]+).*', '\\1', rstr_model@fullmodel),
")) +
ggplot2::geom_bar(stat = 'identity') + ggplot2::ggtitle('",
as.character(get_roi_name(label_desc_df = read_label_desc(rstr_data@roilabeldescfile),roiid=rstr_data@roiids[m])[[1]]),
" ", rstr_data@roimeas,
" vs ", rstr_model@main_effect,"') +
ggplot2::labs(y='",as.character(get_roi_name(label_desc_df = read_label_desc(rstr_data@roilabeldescfile),roiid=rstr_data@roiids[m])[[1]]),
"') +
ggplot2::geom_errorbar(mapping=aes(ymin=len-sd, ymax=len+sd)) +
ggplot2::theme(axis.title=ggplot2::element_text(size=16,face='bold')) +
ggplot2::theme(plot.title=ggplot2::element_text(size=18,face='bold')) +
ggplot2::scale_color_discrete(name= pval_string)\n
ggplot2::ggsave(filename='",
paste0(as.character(get_roi_tag(read_label_desc(rstr_data@roilabeldescfile),rstr_data@roiids[m]))), "_roi",rstr_data@roiids[m],
" ", rstr_data@roimeas,
" vs ", rstr_model@main_effect, ".pdf',device='pdf')\n```\n")
      }
      if (rstr_model@model_type == 'rstr_corr'){
        nb_commands[[m]] <- paste0(sprintf("\n#### Main effect of %s (%d) %s on %s \n",
                                           as.character(get_roi_name(label_desc_df = read_label_desc(rstr_data@roilabeldescfile),roiid=rstr_data@roiids[m])[[1]]),
                                           rstr_data@roiids[m],
                                           rstr_data@roimeas, rstr_model@corr_var),
                                   nb_commands[[m]])
      } else{
        nb_commands[[m]] <- paste0(sprintf("\n#### Main effect of %s (%d) %s on %s controlling for %s \n",
                                           as.character(get_roi_name(label_desc_df = read_label_desc(rstr_data@roilabeldescfile),roiid=rstr_data@roiids[m])[[1]]),
                                           rstr_data@roiids[m],
                                           rstr_data@roimeas, rstr_model@main_effect,rstr_model@covariates ),
                                   nb_commands[[m]])
      }

      }
  } else {
    for (i in 1:length(rstr_data@roiids)){
      current_t_test_table<-paste0("data.frame('P-Values' = ",round(rstr_model@pvalues[i],6),",'T-Values' = ",round(rstr_model@tvalues[i],6),
                                       ",'Adjusted P-Values' = ",round(rstr_model@pvalues_adjusted[i],6), ",'Adjusted T-Values' = ",round(rstr_model@tvalues_adjusted[i],6),
                                   ",row.names = paste0('roiid ',",rstr_data@roiids[i],"))")
      nb_calculations[[i]] <- paste0(sprintf("\n#### T-test output for differences between means of brain imaging phenotypes for %s for roiid %d \n",
                                             rstr_model@group_var,rstr_data@roiids[i]),
                                     nb_calculations[[i]],"\nDT::formatStyle(DT::datatable(",current_t_test_table,
                                     "),column = 'P.Values',color = ifelse(",abs(round(rstr_model@pvalues[i],6)),
                                     "<=0.05,'red','black'))\n```\n")
      nb_plots[[i]]<-paste0(nb_plots[[i]],"ggplot2::ggplot(data=rstr_data@demographics, ggplot2::aes(x=factor(",
                            rstr_model@group_var,"), y = `",selected_col[i],"`, color = factor(",rstr_model@group_var,"))) + ggplot2::geom_boxplot() + ggplot2::ggtitle('",
                            as.character(get_roi_name(label_desc_df = read_label_desc(rstr_data@roilabeldescfile),roiid=rstr_data@roiids[i])[[1]]),
                            " ", rstr_data@roimeas,
                            " vs ", rstr_model@group_var, "') + ggplot2::labs(x = '",rstr_model@group_var,"', y='",
                            as.character(get_roi_name(label_desc_df = read_label_desc(rstr_data@roilabeldescfile),roiid=rstr_data@roiids[i])[[1]]),
                            "', colour = '",rstr_model@group_var,"') +
                            ggplot2::theme(axis.title=ggplot2::element_text(size=16,face='bold')) +
                            ggplot2::theme(plot.title=ggplot2::element_text(size=18,face='bold')) \n
                            ggplot2::ggsave(filename='",
                            paste0(as.character(get_roi_tag(read_label_desc(rstr_data@roilabeldescfile),rstr_data@roiids[i]))), "_roi",rstr_data@roiids[i],
                            "_", rstr_data@roimeas,
                            "_vs_", rstr_model@group_var, ".pdf',device='pdf')\n```\n")
    }
  }



  rmdfileconn<-file(file.path(outdir, sprintf("report_%s_%s.Rmd", rstr_model@model_type, switch(rstr_model@model_type,
                                                                                               rstr_corr = {rstr_model@corr_var},
                                                                                               unpairedttest = {rstr_model@group_var},
                                                                                               pairedttest = {rstr_model@group_var},
                                                                                               rstr_anova = {rstr_model@main_effect},
                                                                                               rstr_lm = {rstr_model@main_effect},
                                                                                               rstr_lmer = {rstr_model@main_effect}))))
  if (rstr_model@model_type == 'rstr_corr'){
    writeLines(c(nb_header, nb_libraries, nb_data_header_one, nb_data_command_one,
                 nb_load_data, nb_data_header_two, nb_data_command_two,
                 nb_data_header_three, nb_data_command_three,
                 unlist(mapply(paste0, nb_commands, nb_plots))), rmdfileconn)
  } else {
    writeLines(c(nb_header, nb_libraries, nb_data_header_one, nb_data_command_one,
                 nb_load_data, nb_data_header_two, nb_data_command_two,
                 nb_data_header_three, nb_data_command_three,
                 unlist(mapply(paste0, nb_commands, nb_calculations, nb_plots))), rmdfileconn)
  }
  close(rmdfileconn)

  # Render the markdown
  rmarkdown::render(file.path(outdir, sprintf("report_%s_%s.Rmd", rstr_model@model_type, switch(rstr_model@model_type,
                                                                                               rstr_corr = {rstr_model@corr_var},
                                                                                               unpairedttest = {rstr_model@group_var},
                                                                                               pairedttest = {rstr_model@group_var},
                                                                                               rstr_anova = {rstr_model@main_effect},
                                                                                               rstr_lm = {rstr_model@main_effect},
                                                                                               rstr_lmer = {rstr_model@main_effect}))),
                    output_file=file.path(outdir, sprintf("report_%s_%s.html", rstr_model@model_type, switch(rstr_model@model_type,
                                                                                                           rstr_corr = {rstr_model@corr_var},
                                                                                                           unpairedttest = {rstr_model@group_var},
                                                                                                           pairedttest = {rstr_model@group_var},
                                                                                                           rstr_anova = {rstr_model@main_effect},
                                                                                                           rstr_lm = {rstr_model@main_effect},
                                                                                                           rstr_lmer = {rstr_model@main_effect}))), quiet = TRUE)

  # Copy modelspec file to the output directory
  file.copy(rstr_model@mspec_file, rstr_out@outdir)
  return(rstr_out)
  }
)

#' Save the statistical analysis output
#' @param rstr_data object of type \code{RstrData}
#' @param rstr_model object of type \code{RstrModel}
#' @param outdir output directory to save the results
#' @param overwrite logical parameter denoting if existing output directory should be overwritten or not (default is FALSE)
#' @param nclusters numeric value denoting number of clusters (default is 10)
#' @export
save_rstr_out <- function(rstr_data, rstr_model, outdir="", overwrite = F, nclusters = 10) {

  outdir <- path.expand(outdir)
  valid_types <- c("sba", "tbm", "roi", "dba", "nca")
  if (! rstr_data@analysis_type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)

  switch(rstr_data@analysis_type,
         sba = { rstr_out <- new("RstrSBAOutput", outdir)},
         tbm = { rstr_out <- new("RstrTBMOutput", outdir) },
         dba = { rstr_out <- new("RstrDBAOutput", outdir) },
         roi = { rstr_out <- new("RstrROIOutput", outdir) }
  )
  if (rstr_data@analysis_type == "tbm" | rstr_data@analysis_type == "dba"){
    rstr_out <- save_out(rstr_out, rstr_data, rstr_model, overwrite = overwrite, nclusters = nclusters)
    invisible(rstr_out)
  # } else if (rstr_data@analysis_type == "sba") {
  #   rstr_out <- save_rstr_out_sba_both_hemi(rstr_out, rstr_data, rstr_model, overwrite = overwrite)
  #   invisible(rstr_out)
  }
  else {
    rstr_out <- save_out(rstr_out, rstr_data, rstr_model, overwrite = overwrite)
    invisible(rstr_out)
  }
}


save_rstr_out_sba_both_hemi <- function(rstr_out, rstr_data, rstr_model, outdir="", overwrite = F) {

  if (rstr_data@hemi == "both") {
    # Split the data array, model and all variables into left and right hemispheres
    save_out(rstr_out, rstr_data, rstr_model, overwrite = overwrite)
    # rstr_data_lh <- rstr_data
    # rstr_data_lh@atlas_filename = rstr_data@atlas_filename_lh
    # rstr_data_lh@atlas_surface <- rstr_data_both@atlas_surface_lh
    # save_out(rstr_out, rstr_data_lh, rstr_model, overwrite = T)
    #
    # rstr_data_rh <- rstr_data
    # rstr_data_rh@atlas_filename = rstr_data@atlas_filename_rh
    # rstr_data_rh@atlas_surface <- rstr_data_both@atlas_surface_rh
    # save_out(rstr_out, rstr_data_rh, rstr_model, overwrite = F)

  }
  else {
    save_out(rstr_out, rstr_data, rstr_model, overwrite = overwrite)
  }
}

#' Save the color LUT and ini colormap images
#' @param measure numeric value denoting the measures used to create the color file
#' @param var_name string denoting name of variable that the color file is being created for
#' @param cmap_title string denoting the type of color map
#' @param rstr_data object of type \code{RstrData}
#' @param rstr_model object of type \code{RstrModel}
#' @param outdir string specifying output directory to save the results in
#' @export

save_rstr_color_files <- function(measure, var_name, cmap_title, rstr_data, rstr_model, outdir) {

  measure <- as.numeric(measure)
  rstr_cmap <- new("RstrColormap", cmap_title, "RdYlBu", measure)
  if (rstr_data@hemi == "both")
    cbar_filename <- paste(paste(rstr_model@model_type, var_name, 'both_hemi', rstr_cmap@cmap_type, sep = '_'), '_cbar.pdf', sep = '')
  else
      cbar_filename <- paste(paste(rstr_model@model_type, var_name, tools::file_path_sans_ext(
      basename(rstr_data@atlas_filename)), rstr_cmap@cmap_type, sep = '_'), '_cbar.pdf', sep = '')
  save_colorbar(file.path(outdir,cbar_filename), rstr_cmap@lut, rstr_cmap@vmin, rstr_cmap@vmax, cmap_title)

  if (rstr_data@hemi == "both")
    cbar_filename <- paste(paste(rstr_model@model_type, var_name, 'both_hemi', rstr_cmap@cmap_type, sep = '_'), '_cbar.png', sep = '')
  else
    cbar_filename <- paste(paste(rstr_model@model_type, var_name, tools::file_path_sans_ext(
      basename(rstr_data@atlas_filename)), rstr_cmap@cmap_type, sep = '_'), '_cbar.png', sep = '')

  save_colorbar(file.path(outdir,cbar_filename), rstr_cmap@lut, rstr_cmap@vmin, rstr_cmap@vmax, cmap_title)

  # save the color LUT
  if (rstr_data@hemi == "both")
    lut_fileprefix <- paste(paste(rstr_model@model_type, var_name, 'both_hemi', rstr_cmap@cmap_type, sep = '_'), '.lut', sep = '')
  else
    lut_fileprefix <- paste(paste(rstr_model@model_type, var_name, tools::file_path_sans_ext(
      basename(rstr_data@atlas_filename)), rstr_cmap@cmap_type, sep = '_'), '.lut', sep = '')

  save_MouseSuiteLUT(file.path(outdir, lut_fileprefix), rstr_cmap@lut)

  # save the cbar json file with colormap for MouseSuite statmap
  #cbarlist <- list("colorbar", rstr_cmap@vmin, rstr_cmap@vmax, t(col2rgb(rstr_cmap@lut)/255))



  cbarlist <- list("colorbar", rstr_cmap@cnegmax, rstr_cmap@cnegmin, rstr_cmap@cposmin, rstr_cmap@cposmax, t(col2rgb(rstr_cmap@lut)/255))
  #names(cbarlist) <- c("jsonid", "cbarmin", "cbarmax", "colormap")
  names(cbarlist) <- c("jsonid", "cbarmin", "cbarlowerthresh", "cbarupperthresh", "cbarmax", "colormap")
  cbar_json <- jsonlite::toJSON(cbarlist, pretty = TRUE, auto_unbox=TRUE)
  if (rstr_data@hemi == "both")
    cbar_json_filename <- paste(paste(rstr_model@model_type, var_name, 'both_hemi', rstr_cmap@cmap_type, sep = '_'), '.cbar', sep = '')
  else
    cbar_json_filename <- paste(paste(rstr_model@model_type, var_name, tools::file_path_sans_ext(
      basename(rstr_data@atlas_filename)), rstr_cmap@cmap_type, sep = '_'), '.cbar', sep = '')

  write(cbar_json, file.path(outdir,cbar_json_filename))

  # save ini colormap with ranges
  if (rstr_data@hemi == "both")
    ini_fileprefix <- paste(paste(rstr_model@model_type, var_name, 'both_hemi', rstr_cmap@cmap_type, sep = '_'), '.ini', sep = '')
  else
    ini_fileprefix <- paste(paste(rstr_model@model_type, var_name, tools::file_path_sans_ext(
      basename(rstr_data@atlas_filename)), rstr_cmap@cmap_type, sep = '_'), '.ini', sep = '')

  save_colormap_to_ini(file.path(outdir, ini_fileprefix), rstr_cmap)
  return(rstr_cmap)
}

#' Save the surface output to the given ouput directory
#' @param measure numeric value denoting the measures used to create the output
#' @param var_name string denoting name of variable used by the function
#' @param rstr_data object of type \code{RstrData}
#' @param rstr_model object of type \code{RstrModel}
#' @param rstr_cmap object of type \code{RstrColormap}
#' @param outdir string specifying output directory to save the results in
#' @export

save_rstr_out_surface <- function(measure, var_name, rstr_cmap, rstr_data, rstr_model, outdir) {

  s1 <- rstr_data@atlas_surface
  s1$attributes <- measure
  s1$vColor <- rstr_cmap@rgbcolors
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  outprefix <- paste(paste(rstr_model@model_type, var_name, tools::file_path_sans_ext(
    basename(rstr_data@atlas_filename)), rstr_cmap@cmap_type, sep = '_'), rstr_data@data_type, sep = '')
  writedfs(file.path(outdir, outprefix), s1)
}

save_rstr_out_surface_both_hemi <- function(measure, var_name, rstr_cmap, rstr_data, rstr_model, outdir) {

  if (rstr_data@hemi == "both") {
    idx_lh <- 1:rstr_data@nvertices_lh
    idx_rh <- (rstr_data@nvertices_lh+1):(rstr_data@nvertices_lh+rstr_data@nvertices_rh)

    rstr_data_lh <- rstr_data
    rstr_data_lh@atlas_filename = rstr_data@atlas_filename_lh
    rstr_data_lh@atlas_surface <- rstr_data@atlas_surface_lh
    measure_lh <- measure[idx_lh]
    rstr_cmap_lh <- rstr_cmap
    rstr_cmap_lh@values <- rstr_cmap@values[idx_lh]
    rstr_cmap_lh@rgbcolors <- rstr_cmap@rgbcolors[idx_lh,1:3]
    save_rstr_out_surface(measure_lh, var_name, rstr_cmap_lh, rstr_data_lh, rstr_model, outdir)

    rstr_data_rh <- rstr_data
    rstr_data_rh@atlas_filename = rstr_data@atlas_filename_rh
    rstr_data_rh@atlas_surface <- rstr_data@atlas_surface_rh
    measure_rh <- measure[idx_rh]
    rstr_cmap_rh <- rstr_cmap
    rstr_cmap_rh@values <- rstr_cmap@values[idx_rh]
    rstr_cmap_rh@rgbcolors <- rstr_cmap@rgbcolors[idx_rh,1:3]
    save_rstr_out_surface(measure_rh, var_name, rstr_cmap_rh, rstr_data_rh, rstr_model, outdir)

  }
  else {
      save_rstr_out_surface(measure, var_name, rstr_cmap, rstr_data, rstr_model, outdir)
  }
}
#' Save the nifti image to the output file
#' @param measure denotes the measure used to create the output
#' @param var_name string denoting name of variable used by the function
#' @param rstr_data object of type \code{RstrData}
#' @param rstr_model object of type \code{RstrModel}
#' @param rstr_cmap object of type \code{RstrColormap}
#' @param outdir string specifying output directory to save the results in
#' @export

save_rstr_out_nifti_image <- function(measure, var_name, rstr_cmap, rstr_data, rstr_model, outdir) {

  outprefix <- paste0(paste(rstr_model@model_type, var_name, tools::file_path_sans_ext(
    basename(rstr_data@atlas_filename)), rstr_cmap@cmap_type, sep = '_'), rstr_data@data_type)
  RNifti::writeNifti(measure, file.path(outdir, outprefix), template = rstr_data@atlas_image)
}

#' Save the measure to the output directory
#' @param measure numeric value denoting the measures used to create the output
#' @param var_name string denoting name of variable used by the function
#' @param label string denoting the label for the object
#' @param rstr_data object of type \code{RstrData}
#' @param rstr_model object of type \code{RstrModel}
#' @param outdir string specifying output directory to save the results in
#' @export

save_rstr_rds <- function(measure, var_name, label, rstr_data, rstr_model, outdir) {

  outprefix <- paste(paste(rstr_model@model_type, var_name, tools::file_path_sans_ext(
    basename(rstr_data@atlas_filename)), label, sep = '_'), ".rds", sep = '')
  saveRDS(measure, file=file.path(outdir, outprefix))
}

#' Save the volume statistics
#' @param log_pvalues log transformed p-values
#' @param log_pvalues_adjusted log transformed adjusted p-values
#' @param tvalues t-values
#' @param tvalues_adjusted adjusted t-values
#' @param var_name string denoting name of variable used by the function
#' @param rstr_data object of type \code{RstrData}
#' @param rstr_model object of type \code{RstrModel}
#' @param outdir string specifying output directory to save the results in
#' @param corr_values correlation values
#' @param corr_values_masked_adjusted adjusted, masked correlation values
#' @export

save_vol_stats_out <- function(log_pvalues, log_pvalues_adjusted, tvalues, tvalues_adjusted,
                                var_name, rstr_data, rstr_model, outdir, corr_values = NULL, corr_values_masked_adjusted = NULL) {

  rstr_cmap <- save_rstr_color_files(log_pvalues, var_name, "log_pvalues", rstr_data, rstr_model, outdir)
  save_rstr_out_nifti_image(log_pvalues, var_name, rstr_cmap, rstr_data, rstr_model, outdir)

  rstr_cmap <- save_rstr_color_files(log_pvalues_adjusted, var_name, "log_pvalues_adjusted", rstr_data, rstr_model, outdir)
  save_rstr_out_nifti_image(log_pvalues_adjusted, var_name, rstr_cmap, rstr_data, rstr_model, outdir)

  rstr_cmap <- save_rstr_color_files(tvalues, var_name, "tvalues", rstr_data, rstr_model, outdir)
  save_rstr_out_nifti_image(tvalues, var_name, rstr_cmap, rstr_data, rstr_model, outdir)

  rstr_cmap <- save_rstr_color_files(tvalues_adjusted, var_name, "tvalues_adjusted", rstr_data, rstr_model, outdir)
  save_rstr_out_nifti_image(tvalues_adjusted, var_name, rstr_cmap, rstr_data, rstr_model, outdir)
  save_rstr_rds(rstr_model@pvalues, var_name, "pvalues", rstr_data, rstr_model, outdir) # Save pvalues as a rds file

  switch(rstr_model@model_type,
         rstr_corr = {
           rstr_cmap <- save_rstr_color_files(corr_values, var_name, "corr_values", rstr_data, rstr_model, outdir)
           save_rstr_out_nifti_image(corr_values, var_name, rstr_cmap, rstr_data, rstr_model, outdir)

           rstr_cmap <- save_rstr_color_files(corr_values_masked_adjusted, var_name, "corr_values_masked_adjusted", rstr_data, rstr_model, outdir)
           save_rstr_out_nifti_image(corr_values_masked_adjusted, var_name, rstr_cmap, rstr_data, rstr_model, outdir)
         }
  )
}

#' Get voxel coordinates of all significant clusters (up to number of clusters)
#' @param rstr_out object of type \code{RstrOut}
#' @param rstr_data object of type \code{RstrData}
#' @param rstr_model object of type \code{RstrModel}
#' @param outdir string specifying output directory to save the results in
#' @param nclusters numeric value specifying number of clusters
#' @export

get_voxelcoord <- function(rstr_out, rstr_data, rstr_model, outdir, nclusters){
  # Call cluster code from terminal
  switch(rstr_model@model_type,
         rstr_anova = {var_name = rstr_model@main_effect},
         rstr_lm = {var_name = rstr_model@main_effect},
         rstr_lmer = {var_name = rstr_model@main_effect},
         rstr_corr = {var_name = rstr_model@corr_var},
         pairedttest = {var_name = rstr_model@group_var},
         unpairedttest = {var_name = rstr_model@group_var}
  )

  if (rstr_model@model_type == "rstr_corr") {
    file_suffix_adjusted <- rs_stat_overlays$corr_values_masked_adjusted
    file_suffix <- rs_stat_overlays$corr_values
  } else {
    file_suffix_adjusted <- rs_stat_overlays$tvalues_adjusted
    file_suffix <- rs_stat_overlays$tvalues
  }

  voxelcoord_call <- paste0("\"",file.path(get_brainsuite_install_path(),rs_binary_files$clustermap),"\" -i \"", outdir,"/",rstr_model@model_type,
                            "_",var_name,"_", tools::file_path_sans_ext(basename(rstr_data@atlas_filename)), "_", file_suffix_adjusted, ".nii.gz\""," -m \"", rstr_data@maskfile, "\" -o \"", outdir, "/cluster.tsv\"", " -n ", nclusters)

  # voxelcoord_call <- paste0("\"",file.path(get_brainsuite_install_path(),rs_binary_files$clustermap),"\" -i \"", outdir,"/",rstr_model@model_type,
  #                           "_",var_name,"_", tools::file_path_sans_ext(basename(rstr_data@atlas_filename)),
  #                           "_tvalues_adjusted.nii.gz\""," -m \"", rstr_data@maskfile, "\" -o \"", outdir, "/cluster.tsv\"", " -n ", nclusters)
  system(voxelcoord_call,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
  if(file.info(paste0(outdir, "/cluster.tsv"))$size!=0){
    vox_table <- read.table(paste0(outdir,"/cluster.tsv"),header=F,sep="\t")
    voxelcoord <- vector("list",nrow(vox_table))
    for (individ_vox in 1:nrow(vox_table)){
      for (vox_component in 1:3){
        voxelcoord[[individ_vox]]<- c(voxelcoord[[individ_vox]],vox_table[individ_vox,3+vox_component])
      }
    }
  }
  # Use tvalues instead of adjusted tvalues if no clusters are found
  if (file.info(paste0(outdir, "/cluster.tsv"))$size==0){
    voxelcoord_call <- paste0("\"",file.path(get_brainsuite_install_path(),rs_binary_files$clustermap),"\" -i \"", outdir,"/",rstr_model@model_type,
                              "_",var_name,"_", tools::file_path_sans_ext(basename(rstr_data@atlas_filename)), "_", file_suffix_adjusted, ".nii.gz\""," -m \"", rstr_data@maskfile, "\" -o \"", outdir, "/cluster.tsv\"", " -n ", nclusters)

#    voxelcoord_call <- paste0("\"",file.path(get_brainsuite_install_path(),rs_binary_files$clustermap),"\" -i \"", outdir,"/",rstr_model@model_type,
#                              "_",var_name,"_", tools::file_path_sans_ext(basename(rstr_data@atlas_filename)), "_tvalues.nii.gz\"",
#                              " -m \"", rstr_data@maskfile, "\" -o \"", outdir, "/cluster.tsv\"", " -n ", nclusters)
    system(voxelcoord_call,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
    if(file.info(paste0(outdir, "/cluster.tsv"))$size!=0){
      vox_table <- read.table(paste0(outdir,"/cluster.tsv"),header=F,sep="\t")
      voxelcoord <- vector("list",nrow(vox_table))
      for (individ_vox in 1:nrow(vox_table)){
        for (vox_component in 1:3){
          voxelcoord[[individ_vox]]<- c(voxelcoord[[individ_vox]],vox_table[individ_vox,3+vox_component])
        }
      }
    } else {voxelcoord <- list(-1)}
  }
  return(voxelcoord)
}

#' Save Rstr Volumetric rmarkdown report including Rmd and the html file
#' This function is used to save Rmd report for TBM
#' @param outdir string specifying output directory to save the results in
#' @param rstr_data object of type `RstrData`
#' @param rstr_model object of type `RstrModel`
#' @param voxelcoord list of peak voxel cluster coordinates
#' @param var_name string specifying the variable name
#' @param stats_string string specifying the statistical measure
#' @export
save_rstr_vol_rmd_html <- function(outdir, rstr_data, rstr_model, voxelcoord, var_name, stats_string){
  
  rmd_preamble <- readLines(system.file("extdata", "report_preamble.Rmd", package="rstr"))
  rmd_preamble <- paste(rmd_preamble, collapse = "\n")
  
  rmd_text <- rmd_preamble
  rstr_report_title_str <- paste0("## ", create_rmd_report_title_str(rstr_data, rstr_model))
  rmd_text <- paste0(rmd_text, rstr_report_title_str, "\n\n")
  # Display the cluster table
  rmd_text <- paste0(rmd_text, "\n\n```{r eval=TRUE, echo=FALSE}\n")
  rmd_text <- paste0(rmd_text, sprintf("vox_table <- read.table(file.path('%s', 'cluster.tsv'), header = F, sep = '	')\n", outdir))
  rmd_text <- paste0(rmd_text, "\ncolnames(vox_table) <- c('Cluster Number', 'Number of Voxels',  'T-Value' , 'X Coord', 'Y Coord', 'Z Coord')\n")
  rmd_text <- paste0(rmd_text, "\nDT::datatable(vox_table, rownames = FALSE)\n")
  rmd_text <- paste0(rmd_text, "```\n\n")
  
  for (jj in seq_along(voxelcoord)) {
    rmd_text <- paste0(rmd_text, sprintf("## Cluster %d: Voxel coordinate (%d, %d, %d)", jj, voxelcoord[[jj]][1], voxelcoord[[jj]][2], voxelcoord[[jj]][3]), " {.tabset}\n")
    for (ii in stats_string){
      rmd_text <- paste0(rmd_text, "\n\n", sprintf("### %s", stats_string_to_rmd_section_title[ii]), "\n")
      rmd_text <- paste0(rmd_text, "\n",  sprintf("```{r cluster%s_%s, fig.cap=''}", jj, ii))
      outprefix <- paste0(rstr_model@model_type, "_", var_name, '_both_hemi_', ii, sep = '_')
      png_filename <- paste0(file.path(outdir,sprintf("cluster%s_", jj)), ii, "_figure.png")
      rmd_text <- paste0(rmd_text, "\n", sprintf("knitr::include_graphics('%s')", png_filename))
      rmd_text <- paste0(rmd_text, "\n```\n")
    }
    rmd_text <- paste0(rmd_text, "\n### {-}\n\n")
  }
  
  rmd_filename <- file.path(outdir, sprintf("report_%s_%s.Rmd", rstr_model@model_type, var_name))
  writeLines(rmd_text, rmd_filename)
  rmarkdown::render(rmd_filename)
  
}