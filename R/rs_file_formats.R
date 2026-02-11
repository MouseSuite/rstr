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

# TODO change the types to adapt to RodentSuite
#' List of file formats used in MouseSuite
#' @export
## TODO: use closures for this in the future
rs_file_formats <- list(
  jacdet = '*.svreg.inv.map.jacdet.*.nii.gz',
  roi_txt = '.roiwise.stats.txt',
  wm_roi_txt = '.wm.stats.tsv',
  svreg_log = '*.svreg.log',
  surf_atlas_left = 'mri.left.mid.cortex.dfs',
  surf_atlas_right = 'mri.right.mid.cortex.dfs',
  nii_atlas = 'mri.bfc.nii.gz',
  nii_maskfile = 'mri.cerebrum.mask.nii.gz'
)

# TODO change the types
#' List ofROI types
#' @export
roi_types <- list(
  svreg_roi_types = c("gmthickness", "gmvolume", "area", "wmvolume"),
  swm_roi_types = c("swmFA", "swmMD", "swmRD", "swmAD"),
  wm_roi_types = c("FA", "FRT_GFA",	"L2",	"L3",	"MD",	"axial",	"mADC",	"radial"),
  gm_roi_types = c("ROI Volume (mm^3)")

)

# TODO change the types to adapt to RodentSuite
#' List of atlas files used in MouseSuite
#' #' @export
#' rs_atlas_files <- list(
#' TODO: Add a list of atlas files supported by MouseSuite
#' )

#' List of binaries used in MouseSuite analysis
#' @export
rs_binary_files <- list(
  clustermap = 'bin/clustermap',
  statmap = 'bin/statmap'
)

#' List of suffixes for atlas files used in MouseSuite
rs_atlas_files_suffix <- list(

  atlas_custom_suffix_tbm = 'bfc.nii.gz',
  atlas_custom_mask_suffix_tbm = 'mask.nii.gz',
  atlas_custom_suffix_dba = 'bfc.nii.gz',
  atlas_custom_mask_suffix_dba = 'wm.mask.nii.gz'
)

analysis_type_list <- list(
  sba = 'sba',
  tbm = 'tbm',
  roi = 'roi',
  dba = 'dba',
  nca = 'nca'
)

#' Mapping of analysis_type to readable string
rs_analysis_type_to_readable_string <- list(
  tbm = "TBM",
  roi = "ROI",
  dba = "DBA"
)

#' Mapping of stats measures to readable string to be used as a section title on the Rmd report
stats_string_to_rmd_section_title <- list(
  log_pvalues_adjusted = "Adjusted P-values",
  log_pvalues = "P-values",
  tvalues_adjusted = "Adjusted T-values",
  tvalues = "T-values",
  corr_values_masked_adjusted = "Masked Adjusted correlations",
  corr_values = "Correlations"
)

rs_data_types <- list(
  surface = '.dfs',
  nifti_image = '.nii.gz'
)

#' List of statistical and other overlay types used in MouseSuite
rs_stat_overlays <- list(
  log_pvalues_adjusted = "log_pvalues_adjusted",
  tvalues_adjusted = "tvalues_adjusted",
  log_pvalues = "log_pvalues",
  tvalues = "tvalues",
  pvalues = "pvalues",
  corr_values = "corr_values",
  corr_values_masked_adjusted = "corr_values_masked_adjusted"
)
#' Returns the filename for the designated image
#' @param outdir string denoting the output directory
#' @param voxelcoord all outputted voxelcoordinates
#' @param overlay_name string denoting the name of the overlay
#' @param brain_sector_index numeric value denoting which brain slice is desired
#' @param voxelcoord_index numeric value denoting desired voxel coordinate
#'
get_render_image_filename <- function(outdir,voxelcoord, overlay_name, brain_sector_index, voxelcoord_index) {
  view_order <- c("sag","cor","ax")
  return(paste0(outdir, "/png_images_crosshairs/", view_order[brain_sector_index], voxelcoord[[voxelcoord_index]][brain_sector_index],"_",overlay_name,"_cluster",voxelcoord_index,".png"))
}
#' Generate the designated cortical surface file
#' @param hemi string denoting which hemisphere of the brain is desired
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
rs_surface_file_string <- function(hemi="left", smooth = 0) {

  if (smooth != 0)
    return(paste('atlas.pvc-thickness_0-6mm.', sprintf('smooth%2.1fmm.', smooth), hemi, '.mid.cortex.dfs', sep = ''))
  else
    return(paste('atlas.pvc-thickness_0-6mm.', hemi, '.mid.cortex.dfs', sep = ''))
}
#' Generate the designated tensor-based morphometry file
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
rs_volume_jacobian_file_string <- function(smooth = 0) {

  if (smooth != 0)
    return(paste('%s.svreg.inv.jacobian.', sprintf('smooth%2.1fmm.nii.gz', smooth), sep = ''))
  else
    return('%s.svreg.inv.jacobian.nii.gz')
}

#' Generate the designated BIDS compatible (T1w) tensor-based morphometry file
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
rs_BIDS_volume_jacobian_file_string <- function(smooth = 0) {

  if (smooth != 0)
    return(paste('%s_T2w.gsmooth.inv.warp-Jacobian-', sprintf('smooth%2.1fmm.nii.gz', smooth), sep = ''))
  else
    return("%s_T2w.inv.warp-Jacobian.nii.gz")
}

#' Generate the designated diffusion surface file
#' @param measure string designating the type of diffusion measure used
#' @param smooth numeric value designating the smoothing used (default is 0)
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#'
rs_diffusion_file_string <- function(measure = "FA", smooth = 0, eddy = TRUE) {

  valid_diffusion_measures <- c('FA', 'MD', 'axial', 'radial', 'mADC', 'FRT_GFA')
  if (! measure %in% valid_diffusion_measures)
    stop(sprintf('Invalid diffusion measure: %s. Valid measures are %s.', measure, paste(valid_diffusion_measures, collapse = ', ')))

  if (eddy == TRUE && smooth != 0)
    return(paste0('%s.dwi.RAS.correct.atlas.', measure, sprintf('.smooth%2.1fmm.nii.gz', smooth)))
  if (eddy == FALSE && smooth != 0)
    return(paste0('%s.dwi.RAS.atlas.', measure, sprintf('.smooth%2.1fmm.nii.gz', smooth)))
  if (eddy == TRUE && smooth == 0)
    return(paste0('%s.dwi.RAS.correct.atlas.', measure, '.nii.gz'))
  if (eddy == FALSE && smooth == 0)
    return(paste0('%s.dwi.RAS.atlas.', measure, '.nii.gz'))
}

#' Generate the BIDS designated diffusion surface file
#' @param measure string designating the type of diffusion measure used
#' @param smooth numeric value designating the smoothing used (default is 0)
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#'
rs_BIDS_diffusion_file_string <- function(measure = "FA", smooth = 0, eddy = TRUE) {

  valid_diffusion_measures <- c('FA', 'MD', 'axial', 'radial', 'mADC', 'FRT_GFA')
  if (! measure %in% valid_diffusion_measures)
    stop(sprintf('Invalid diffusion measure: %s. Valid measures are %s.', measure, paste(valid_diffusion_measures, collapse = ', ')))

  if (eddy == TRUE && smooth != 0)
    return(paste0('%s_dwi.dwi.RAS.correct.atlas.', measure, sprintf('.smooth%2.1fmm.nii.gz', smooth)))
  if (eddy == FALSE && smooth != 0)
    return(paste0('%s_dwi.dwi.RAS.atlas.', measure, sprintf('.smooth%2.1fmm.nii.gz', smooth)))
  if (eddy == TRUE && smooth == 0)
    return(paste0('%s_dwi.dwi.RAS.correct.atlas.', measure, '.nii.gz'))
  if (eddy == FALSE && smooth == 0)
    return(paste0('%s_dwi.dwi.RAS.atlas.', measure, '.nii.gz'))
}

#' Stops analysis if desired type of analysis is not a valid analysis type
#' @param analysis_type string denoting desired type of analysis to be performed
#'
get_rs_file_list <- function(analysis_type) {
  valid_analysis_types <- unlist(analysis_type_list, use.names = FALSE)
  if (!(analysis_type %in% valid_analysis_types)) {
    stop(sprintf('Valid brain analyses are %s', paste(unlist(analysis_type_list), collapse = ', ')),
         call. = FALSE)
  }


}
#' Returns a list of all ROI files for all subjects
#' @param rstr_data object of type \code{RstrData}
#' @param roistat_fileprefix file prefix for the ROI stats file (either .tissue.stats.txt or .roiwise.stats.txt)",
#' @param roilabeldescfile an xml file containing ROI label descriptions",
#'
get_roi_file_list <- function(rstr_data, roistat_fileprefix, roilabeldescfile) {

  bids_flag_and_filelist <- check_bids_compatibility_and_get_filelist(rstr_data, type="roi", hemi="", smooth="", measure="", roistat_fileprefix=roistat_fileprefix, roilabeldescfile=roilabeldescfile)
  return(bids_flag_and_filelist$filelist)

}
#' Returns a list of the cortical surface files for all subjects
#' @param rstr_data object of type \code{RstrData}
#' @param hemi designates which hemisphere is of interest
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
get_sba_file_list <- function(rstr_data, hemi, smooth = 0) {

  bids_flag_and_filelist <- check_bids_compatibility_and_get_filelist(rstr_data, type="sba", hemi, smooth)
  return(bids_flag_and_filelist$filelist)
}
#' Returns a list of the tensor-based files for all subjects
#' @param rstr_data object of type \code{RstrData}
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
get_tbm_file_list <- function(rstr_data, smooth = 0) {

  bids_flag_and_filelist <- check_bids_compatibility_and_get_filelist(rstr_data, type="tbm", hemi="left", smooth)
  return(bids_flag_and_filelist$filelist)
}

#' Returns a list of the diffusion files for all subjects
#' @param rstr_data object of type \code{RstrData}
#' @param measure numeric value denoting the measures used to create the output
#' @param smooth numeric value designating the smoothing used (default is 0)
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#'
get_dba_file_list <- function(rstr_data, measure, smooth = 0, eddy = TRUE) {

  bids_flag_and_filelist <- check_bids_compatibility_and_get_filelist(rstr_data, type="dba", hemi="left", smooth, measure)
  return(bids_flag_and_filelist$filelist)
}
#' Read the MouseSuite atlas prefix from the atlas
#' @param atlas filepath for atlas
#'
get_mousesuite_custom_volume_atlas_prefix <- function(atlas) {

  # Check if the parent directory exists for the atlas
  check_file_exists(dirname(atlas), raise_error = TRUE)

  # If atlas points to a nifti image, return the prefix of the atlas (everything until the bfc.nii.gz)
  if (substr(atlas, nchar(atlas) - 9, nchar(atlas)) == rs_atlas_files_suffix$atlas_custom_suffix_tbm)
    return(substr(atlas, 1, nchar(atlas) - 11))

  # If atlas points to a valid prefix, return atlas
  if (!identical(Sys.glob(file.path(paste(atlas, "*", rs_atlas_files_suffix$atlas_custom_suffix_tbm, sep = ""))), character(0)))
    return(atlas)

  # Otherwise raise an exception
  stop(sprintf('Invalid custom atlas prefix/path: %s', atlas), call. = FALSE)

}
#' Read the MouseSuite atlas path from the svreg log file.
#' @param logfile path to the log file present in an individual subject directory
#'
get_mousesuite_atlas_path_from_logfile <- function(logfile) {
# TODO: In the future, support loading the atlas from a log file
}

#' Read the MouseSuite atlas identifier from the svreg log file.
#' @param logfile path to the log file present in an individual subject directory
#'
get_mousesuite_atlas_id_from_logfile <- function(logfile) {
  fid = file(logfile, "rt")
  log_lines <- readLines(fid, n=2)
  close(fid)
  # TODO: Add the functionality for loading the log file
}

#' Gets the MouseSuite log file's filepath
#' @param subjdir individual subject directory that the log file exists in
#' @param csv csv file for the log file
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @export
get_mousesuite_logfilename <- function(subjdir, csv, exclude_col="") {
  # Open the log file and get the atlas file name
  if ( identical(tools::file_ext(csv), 'tsv') | identical(tools::file_ext(csv), 'csv')  ) {
    demog <- read_demographics(csv, exclude_col)
  }
  # The first column has to contain subject IDs which are same as subject directories
  first_subjid <- demog[[1]][1]
  # Get atlas names from log files.
  log_file <- Sys.glob(file.path(subjdir, first_subjid, "*", '*.log')) #Extra * due to BIDS compatibility (extra dir anat or dwi)

  if (length(log_file) == 0){
    log_file <- Sys.glob(file.path(subjdir, first_subjid, '*.log'))  #If subject data is not BIDS compatible
  }
  # log_file <- file.path(subjdir, first_subjid, sprintf('%s.log', first_subjid))
  if (check_file_exists(log_file, raise_error = TRUE,
                    errmesg = sprintf('Could not find log in the subject directory %s/%s. Please check if the subject directory is valid.', subjdir, first_subjid)))
    return(tools::file_path_as_absolute(log_file))
}

#' Get the svreg log file for each subject
#' @param subjdir individual subject directory that the svreg.log file exists in
#' @param csv csv file for the svreg.log file
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @export
get_mousesuite_logfilename_for_all_subjects <- function(subjdir, csv, exclude_col="") {
  # Open the svreg.log file and get the atlas file name
  if ( identical(tools::file_ext(csv), 'tsv') | identical(tools::file_ext(csv), 'csv')  ) {
    demog <- read_demographics(csv, exclude_col)
  }
  # The first column has to contain subject IDs which are same as subject directories
  first_subjid <- demog[[1]][1]
  # Get atlas names from log files.
  log_files <- Sys.glob(file.path(subjdir, demog[[1]], "*", '*.svreg.log'))

  if (length(log_files) == 0){
    log_files <- Sys.glob(file.path(subjdir, demog[[1]], '*.svreg.log'))
  }
  # log_file <- file.path(subjdir, first_subjid, sprintf('%s.svreg.log', first_subjid))
  if (check_multiple_files_exists(log_files, errmesg = 'Could not find svreg.log in a few subject directories.')) {
    return(log_files)
  }
}

#' Get the desired cortical surface atlas
#' @param mousesuite_atlas_id individual subject directory that the log file exists in
#' @param hemi designates which hemisphere of the brain
#'
get_sba_atlas <- function(mousesuite_atlas_id, hemi) {

# TODO: In the future support cortical surface based atlases in MouseSuite
}

#' Get the tensor based morphometry atlas and mask
#' @param mousesuite_atlas_id individual subject directory that the svreg.log file exists in
#'
get_tbm_atlas_and_mask <- function(mousesuite_atlas_id) {

# TODO: In the future support multiple volumetric atlases in MouseSuite
  
}

#' Check and return the maskfile and the atlas file
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param maskfile filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
check_tbm_atlas_and_mask <- function(subjdir="", csv="", atlas="", maskfile="", exclude_col) {

  tbm_atlas_and_mask <- list("nii_atlas" = atlas, "nii_atlas_mask" = maskfile)
  if (maskfile == "")
    tbm_atlas_and_mask$nii_atlas_mask <- maskfile

  check_file_exists(atlas, raise_error = TRUE)
  tbm_atlas_and_mask <- list("nii_atlas" = atlas, "nii_atlas_mask" = maskfile)
  return(tbm_atlas_and_mask)
  
  
  # if (maskfile == ""  && atlas == "") {
  #   mousesuite_atlas_id <- get_mousesuite_atlas_id_from_logfile(get_mousesuite_logfilename(subjdir, csv, exclude_col))
  #   tbm_atlas_and_mask <- get_tbm_atlas_and_mask(mousesuite_atlas_id)
  # }
  # else {
  #   if (maskfile != "" && atlas == "") {
  #     mousesuite_atlas_id <- get_mousesuite_atlas_id_from_logfile(get_mousesuite_logfilename(subjdir, csv, exclude_col))
  #     check_file_exists(maskfile, raise_error = TRUE)
  #     tbm_atlas_and_mask$nii_atlas_mask <- maskfile
  #     tbm_atlas_and_mask$nii_atlas <- get_tbm_atlas(mousesuite_atlas_id)
  #   }
  #   else if (maskfile== "" && atlas != "") {
  #     mousesuite_atlas_id <- get_mousesuite_atlas_id_from_logfile(get_mousesuite_logfilename(subjdir, csv, exclude_col))
  #     check_file_exists(atlas, raise_error = TRUE)
  #     tbm_atlas_and_mask$nii_atlas <- atlas
  #     tbm_atlas_and_mask$nii_atlas_mask <- get_tbm_mask(mousesuite_atlas_id)
  #   }
  #   else if (maskfile!= "" && atlas != "") {
  #     check_file_exists(maskfile, raise_error = TRUE)
  #     check_file_exists(atlas, raise_error = TRUE)
  #     tbm_atlas_and_mask <- list("nii_atlas" = atlas, "nii_atlas_mask" = maskfile)
  #   }
  # }
  # 
  # return(tbm_atlas_and_mask)

}



#' Get the MouseSuite tensor based morphometry atlas
#' @param mousesuite_atlas_id individual subject directory that the log file exists in
#'
get_tbm_atlas <- function(mousesuite_atlas_id) {
# TODO: Support loading volumetric atlas

}

#' Get the MouseSuite tensor based morphometry mask
#' @param mousesuite_atlas_id individual subject directory that the svreg.log file exists in
#'
get_tbm_mask <- function(mousesuite_atlas_id) {
# TODO: Support loading volumetric masks
}


#' Check that cortical surface atlas exists
#' @param atlas filepath for sba atlas
#' @param maskfile filepath for the atlas mask file
#'
get_custom_sba_atlas_and_mask <- function(atlas, maskfile="") {
  #TODO Only the atlas file is implemented
  check_file_exists(atlas, raise_error = TRUE)
  return(atlas)
}
#' Get the custom tensor based morphometry atlas and mask
#' @param mousesuite_custom_atlas_prefix file prefix for atlas
#'
get_custom_tbm_atlas_and_mask <- function(mousesuite_custom_atlas_prefix) {

  mousesuite_custom_atlas_prefix <- get_mousesuite_custom_volume_atlas_prefix(mousesuite_custom_atlas_prefix)
  nii_atlas <- paste(mousesuite_custom_atlas_prefix, ".", rs_atlas_files_suffix$atlas_custom_suffix_tbm, sep="")
  nii_atlas_mask <- paste(mousesuite_custom_atlas_prefix, ".", rs_atlas_files_suffix$atlas_custom_mask_suffix_tbm, sep="")
  check_file_exists(nii_atlas, raise_error = TRUE)
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(list("nii_atlas" = nii_atlas, "nii_atlas_mask" = nii_atlas_mask))
}
#' Get the diffusion atlas and mask
#' @param mousesuite_atlas_id individual subject directory that the svreg.log file exists in
#'
get_dba_atlas_and_mask <- function(mousesuite_atlas_id) {
# TODO: In the future, support loading rodent diffusion atlas and mask
}

#' Get the MouseSuite tensor based morphometry atlas
#' @param mousesuite_atlas_id individual subject directory that the svreg.log file exists in
#'
get_dba_atlas <- function(mousesuite_atlas_id) {
# TODO: In the future, support loading rodent diffusion atlas
}

#' Get the MouseSuite tensor based morphometry mask
#' @param mousesuite_atlas_id individual subject directory that the svreg.log file exists in
#'
get_dba_mask <- function(mousesuite_atlas_id) {
# TODO: In the future, support loading rodent diffusion mask
}



#' Get the diffusion custom atlas and mask
#' @param mousesuite_custom_atlas_prefix file prefix for atlas
#'
get_custom_dba_atlas_and_mask <- function(mousesuite_custom_atlas_prefix) {

  mousesuite_custom_atlas_prefix <- get_mousesuite_custom_volume_atlas_prefix(mousesuite_custom_atlas_prefix)
  nii_atlas <- paste(mousesuite_custom_atlas_prefix, ".", rs_atlas_files_suffix$atlas_custom_suffix_dba, sep = "")
  nii_atlas_mask <- paste(mousesuite_custom_atlas_prefix, ".", rs_atlas_files_suffix$atlas_custom_mask_suffix_dba, sep = "")
  check_file_exists(nii_atlas, raise_error = TRUE)
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(list("nii_atlas" = nii_atlas, "nii_atlas_mask" = nii_atlas_mask))
}
#' Reads the demographics from an inputted csv file
#' @param csvfile csv file containing the demographics
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @export
read_demographics <- function(csvfile, exclude_col="") {

  switch(tools::file_ext(csvfile),
         "tsv" = {demog <- read.table(file = csvfile, sep = "\t", header = T)},
         "csv" = {demog <- read.csv(csvfile)})
  # demo <- read.csv(csvfile)
  colnames(demog)[1] <- "subjID"
  if (exclude_col != "") {
    if (! exclude_col %in% colnames(demog))
      stop(sprintf("Exclude column specified as %s does not exist in %s.", exclude_col, csvfile), call. = FALSE)
    demog <- demog[demog[, exclude_col] == 0,]
  }
  return(demog)
}

check_bids_compatibility_and_get_filelist <- function(rstr_data, type="sba", hemi, smooth = 0, measure="FA", eddy = TRUE, roimeas = "ROI Volume (mm^3)",
                                                      roistat_fileprefix, roilabeldescfile) {

  valid_types <- c("sba", "tbm", "roi", "dba")
  if (! type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)

  if (type == "sba") {
    filelist <- file.path(rstr_data@subjdir, rstr_data@demographics$subjID, rs_surface_file_string(hemi, smooth))
    bids_filelist <- file.path(rstr_data@subjdir, rstr_data@demographics$subjID, 'anat', rs_surface_file_string(hemi, smooth))
  }
  else if (type == "tbm") {
    filelist <- file.path(rstr_data@subjdir, rstr_data@demographics$subjID, sprintf(rs_volume_jacobian_file_string(smooth), rstr_data@demographics$subjID))
    bids_filelist <- file.path(rstr_data@subjdir, rstr_data@demographics$subjID, 'anat', sprintf(rs_BIDS_volume_jacobian_file_string(smooth), rstr_data@demographics$subjID))
  }
  else if (type == "dba") {
    valid_dba_measures <- c('FA', 'MD', 'axial', 'radial', 'mADC', 'FRT_GFA')
    if (! measure %in% valid_dba_measures) {
      stop(sprintf("Valid dba measures are %s.", paste(valid_dba_measures, collapse = ', ')), call. = FALSE)
    }
    filelist <- file.path(rstr_data@subjdir, rstr_data@demographics$subjID, sprintf(rs_diffusion_file_string(measure, smooth, eddy), rstr_data@demographics$subjID))
    bids_filelist <- file.path(rstr_data@subjdir, rstr_data@demographics$subjID, 'dwi', sprintf(rs_BIDS_diffusion_file_string(measure, smooth, eddy), rstr_data@demographics$subjID))

  }
  else if (type == "roi") {
    valid_roi_measures <- c(roi_types$svreg_roi_types, roi_types$swm_roi_types, roi_types$swm_roi_types, roi_types$gm_roi_types)
    if (roimeas %in% roi_types$gm_roi_types) {
      bids_filelist <- file.path(rstr_data@subjdir, rstr_data@demographics$subjID, 'anat',
                                     sprintf('%s%s%s', rstr_data@demographics$subjID, '_T2w', roistat_fileprefix))
      filelist <- file.path(rstr_data@subjdir, rstr_data@demographics$subjID,
                                 sprintf('%s%s%s', rstr_data@demographics$subjID, '_T2w', roistat_fileprefix))

    }
    else if (roimeas %in% roi_types$wm_roi_types) {
      bids_filelist <- file.path(rstr_data@subjdir, rstr_data@demographics$subjID, 'dwi',
                                     sprintf('%s%s%s', rstr_data@demographics$subjID, '_T2w', roistat_fileprefix))
      filelist <- file.path(rstr_data@subjdir, rstr_data@demographics$subjID,
                                 sprintf('%s%s%s', rstr_data@demographics$subjID, '_T2w', roistat_fileprefix))

    }
  }

  # Check BIDS compatibility
  if (any(file.exists(filelist))) {
    bids_compatible = FALSE
  }
  else if (any(file.exists(bids_filelist))) {
    bids_compatible = TRUE
    filelist <- bids_filelist
  }
  else {
    message('Following subjects have missing files')
    print(filelist[which(!file.exists(filelist))], row.names = FALSE)
    stop('\nCheck if SVREG was run succesfully on all the subjects. Optionally, check the smoothing level (smooth= under [subject]).\nIt is possible that surface or volume files at the specified smoothing level do not exist.',
         call. = FALSE)

#    stop('\nCould not understand the directory hierarchy. Check if svreg was run succesfully on all the subjects. Also check the smoothing level (smooth= under [subject]).\nIt is possible that surface files at the specified smoothing level do not exist.',
#         call. = FALSE)
  }

  # Check if all files exist
  if ( !all(file.exists(filelist))) {
    message('Following subjects have missing files')
    print(filelist[which(!file.exists(filelist))], row.names = FALSE)
    stop('\nCheck if SVREG was run succesfully on all the subjects. Optionally, check the smoothing level (smooth= under [subject]).\nIt is possible that surface or volume files at the specified smoothing level do not exist.',
         call. = FALSE)
  }
  return (list("bids_compatible" = bids_compatible, "filelist" = filelist))

}


create_rmd_report_title_str <- function(rstr_data, rstr_model) {
  
  if(rstr_data@analysis_type == "dba") {
    title_string <- paste(rstr:::rs_analysis_type_to_readable_string[[rstr_data@analysis_type]], " ", rstr_data@measure, " ",
                          rstr:::rs_model_type_to_readable_string[[rstr_model@model_type]], '-- ', sep="")
  }
  else {
    title_string <- paste(rstr:::rs_analysis_type_to_readable_string[[rstr_data@analysis_type]], " ",
                          rstr:::rs_model_type_to_readable_string[[rstr_model@model_type]], '-- ', sep="")
  }
  
  switch(rstr_model@model_type,
         rstr_anova = {
           title_string <- paste(title_string, "Main effect of", rstr_model@main_effect, "with covariates",
                                 rstr_model@covariates)
         },
         rstr_lm = {
           title_string <- paste(title_string, "Main effect of", rstr_model@main_effect, "with covariates",
                                 rstr_model@covariates)
         },
         rstr_lmer = {
           title_string <- paste(title_string, "Fixed effect of", rstr_model@main_effect, "with covariates",
                                 rstr_model@covariates, ", Random effect of", rstr_model@group_var)
         },
         rstr_corr = {
           title_string <- paste(title_string, "with", rstr_model@corr_var)
         },
         pairedttest = {
           title_string <- paste(title_string, "between samples for", rstr_model@group_var)
         },
         unpairedttest = {
           title_string <- paste(title_string, "between samples for", rstr_model@group_var)
         }
  )
  
  return(title_string)
}