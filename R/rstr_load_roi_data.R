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

#' Load ROI data from subjects (BIDS compatible)
#'
#' Loads the ROI measure for statistical analysis.
#' @param subjects_dir character string for subject directory.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' @param roiids vector of numeric label identifiers for the region of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. 
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @param roistat_fileprefix file prefix for the ROI stats file (either .tissue.stats.txt or .roiwise.stats.txt)",
#' @param roilabeldescfile an xml file containing ROI label descriptions",
#' @export
#'
rstr_load_roi_data <- function(subjects_dir, csv, roiids, roimeas, exclude_col, roistat_fileprefix, roilabeldescfile) {

  if (!dir.exists(subjects_dir)) {
    stop(sprintf("Subjects directory %s does not exist.\n", subjects_dir), call. = FALSE)
  }

  demographics <- read_demographics(csv, exclude_col)
  if("File_roi" %in% colnames(demographics)){
    warning(sprintf("The file %s already contains a File_roi column.\nWill overwrite this column.\n", csv), call.=FALSE)
  }

  demographics$subjID <- as.character(demographics$subjID)

  rstr_data <- new("RstrROIData", subjects_dir, csv, exclude_col)
  roiwise_file_list <- get_roi_file_list(rstr_data, roistat_fileprefix, roilabeldescfile)

  if (any(is.na(roiwise_file_list))) {
    stop(sprintf("Subject %s is either missing the roiwise.stats.txt or the wm.stats.tsv file depending on the ROI measure used.\n",
                 demographics$subjID[which(is.na(roiwise_file_list))]), call. = FALSE)
  }

  demographics$File_roi <- unlist(roiwise_file_list)

  roi_data_frame <- read_roi_data_for_all_subjects(demographics$File_roi, demographics,roiids, roimeas, roistat_fileprefix, roilabeldescfile)

  # Put outputted data frame together with demographics data frame
  combined_roidata_and_demographics <- cbind(demographics[,-which(names(demographics) == "File_roi")], roi_data_frame)

  rstr_data <- list(df=combined_roidata_and_demographics, roiids=roiids, roimeas=roimeas)

  # Finally also include the command to load the data
  pasted_roiids <- paste(roiids,collapse = ", ")
  rstr_data$load_data_command <- sprintf("rstr_data <- rstr_load_roi_data( '%s', '%s', c( %s), '%s', '%s') ",
                                        subjects_dir, csv, pasted_roiids, roimeas, exclude_col)

  return(rstr_data)
}



