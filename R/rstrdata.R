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

#' S4 class for storing data for statistical analysis
#' @slot data_array matrix containing data of dimensions (N x T), where N = number of subjects and T = number of vertices/voxels.
#' @slot data_array_lh matrix containing data for left hemisphere.
#' @slot data_array_rh matrix containing data for right hemisphere.
#' @slot analysis_type character string denoting the type of analysis. Valid types are "sba", "tbm" or "roi".
#' @slot data_type character string denoting the type of data. Valid types are "surface" or "nifti_image".
#' @slot demographics data.frame containing the demographic information. Usually loaded from a csv file.
#' @slot subjdir character string for subject directory.
#' @slot csv filename of a comma separated (csv) file containing the subject demographic information.
#' @slot smooth numeric value used to smooth the data.
#' @slot measure character string denoting the type of measure used.
#' @slot filelist list of files belonging to N subjects.
#' @slot load_data_command character string for the command used to load the data
#' @slot hemi character string to specifiy the hemisphere. Valid values are "left" or "right" or "both"
#' @slot nvertices_lh numeric value containing the number of vertices for the left cortical surface atlas
#' @slot nvertices_rh numeric value containing the number of vertices for the right cortical surface atlas
#'
#' @export
RstrData <- setClass(
  "RstrData",
  slots = list(
    data_array = "matrix",
    data_array_lh = "matrix",
    data_array_rh = "matrix",
    analysis_type = "character",
    data_type = "character",
    demographics = "data.frame",
    subjdir = "character",
    csv = "character",
    smooth = "numeric",
    measure = "character",
    filelist = "character",
    load_data_command = "character",
    hemi = "character",
    nvertices_lh = "numeric",
    nvertices_rh = "numeric"
  )
)

RstrROIData <- setClass(
  "RstrROIData",
  slots = list(roiids = "numeric",
               roimeas = "character",
               roilabeldescfile = "character",
               roistat_fileprefix = "character"),
  contains = "RstrData"
)

RstrSBAData <- setClass(
  "RstrSBAData",
  slots = list(atlas_filename = "character",
               atlas_filename_lh = "character",
               atlas_filename_rh = "character",
               atlas_surface = 'list',
               atlas_surface_lh = 'list',
               atlas_surface_rh = 'list'),
  contains = "RstrData"
)

setOldClass("niftiImage") #Declare niftiImage (RNifti) so the @slot atlas_image can be defined

RstrTBMData <- setClass(
  "RstrTBMData",
  slots = list(atlas_filename = "character",
               atlas_image = 'niftiImage',
               maskfile = 'character',
               mask_idx = 'vector'),
  contains = "RstrData"
)

RstrDBAData <- setClass(
  "RstrDBAData",
  slots = list(atlas_filename = "character",
               atlas_image = 'niftiImage',
               maskfile = 'character',
               mask_idx = 'vector'),
  contains = "RstrData"
)

setMethod("initialize", valueClass = "RstrData", signature = "RstrData", function(.Object, subjdir, csv, exclude_col) {
  
  check_file_exists(subjdir, raise_error = TRUE)
  check_file_exists(csv, raise_error = TRUE)
  .Object@subjdir = subjdir
  .Object@csv <- csv
  .Object@demographics <- read_demographics(csv, exclude_col = exclude_col)
  return(.Object)
})

#' A generic function to load data for statistical analysis.
#' @param rstr_data object of type \code{\link{RstrData}}
#' @param atlas_filename path name to the atlas
#' @param maskfile path name to the mask file
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right" or "both".
#' @param measure character specifying the brain imaging measure. If analyzing diffusion data, should be "FA".
#' @param smooth numeric value denoting the smoothing level.
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#' @param roiids numeric label identifier for the region of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. 
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @param roistat_fileprefix file prefix for the ROI stats file (either .tissue.stats.txt or .roiwise.stats.txt)",
#' @param roilabeldescfile an xml file containing ROI label descriptions",
#' @details
#' For the most part, the user will never have to call this function directly.
#' Instead the user should call \code{\link{load_rstr_data}}.
#' @seealso \code{\link{load_rstr_data}}
#'
#' @export
setGeneric("load_data", valueClass = "RstrData", function(rstr_data, atlas_filename = NULL, maskfile = NULL, hemi = "left", measure = "", smooth = 0.0, eddy = TRUE, roiids = NULL, roimeas = NULL, exclude_col, roistat_fileprefix, roilabeldescfile) {
  standardGeneric("load_data")
})

setGeneric("load_demographics", valueClass = "RstrData", function(object) {
  standardGeneric("load_demographics")
})

#' @rdname load_data
setMethod("load_data", signature = "RstrData", function(rstr_data, roiids = NULL, roimeas = NULL, exclude_col, roistat_fileprefix, roilabeldescfile) {
  return(rstr_data)
})

#' @rdname load_data
setMethod("load_data", signature = "RstrSBAData", function(rstr_data, atlas_filename, hemi, smooth) {
  
  rstr_data@atlas_filename <- atlas_filename
  rstr_data@atlas_surface <- readdfs(atlas_filename)
  sba_filelist <- get_sba_file_list(rstr_data, hemi, smooth)
  attrib_siz <- rstr_data@atlas_surface$hdr$nVertices
  rstr_data@data_array <- read_dfs_attributes_for_all_subjects(sba_filelist, attrib_siz)
  rstr_data@filelist <- sba_filelist
  rstr_data@analysis_type <- "sba"
  rstr_data@smooth <- smooth
  rstr_data@data_type <- rs_data_types$surface
  return(rstr_data)
})

#' @rdname load_data
setMethod("load_data", signature = "RstrTBMData", function(rstr_data, atlas_filename, maskfile = NULL, smooth) {
  
  rstr_data@atlas_filename <- atlas_filename
  rstr_data@atlas_image <- RNifti::readNifti(atlas_filename)
  rstr_data@filelist <- get_tbm_file_list(rstr_data, smooth)
  rstr_data@smooth <- smooth
  rstr_data@measure <- ""
  attrib_siz <- length(rstr_data@atlas_image)
  if ( !is.null(maskfile) && maskfile != "" ) {
    rstr_data@maskfile <- maskfile
    mask_image <- as.vector(RNifti::readNifti(maskfile))
    if ( length(mask_image) != attrib_siz) {
      stop(sprintf('Dimensions of atlas file %s and maskfile %s do not match', atlas_filename, maskfile), call. = FALSE)
    }
    rstr_data@mask_idx <- which(mask_image > 0)
  }
  else
    rstr_data@mask_idx = 1:attrib_siz
  
  first_subject_file <- as.vector(RNifti::readNifti(rstr_data@filelist[1]))
  # Check if the dimensions of first subject file and atlas match (the dimensions of atlas and mask are already checked above)
  if ( length(first_subject_file) != attrib_siz) {
    stop(sprintf('Dimensions of the atlas file %s and subject %s do not match. Check if you are using the correct atlas', atlas_filename, rstr_data@filelist[1]), call. = FALSE)
  }
  
  rstr_data@data_array <- read_nii_images_for_all_subjects(rstr_data@filelist, attrib_siz, rstr_data@mask_idx)
  rstr_data@analysis_type <- "tbm"
  rstr_data@data_type <- rs_data_types$nifti_image
  rstr_data@hemi <- "NA"
  return(rstr_data)
})

#' @rdname load_data
setMethod("load_data", signature = "RstrDBAData", function(rstr_data, atlas_filename, maskfile = NULL, measure, smooth, eddy) {
  
  rstr_data@atlas_filename <- atlas_filename
  rstr_data@atlas_image <- RNifti::readNifti(atlas_filename)
  rstr_data@filelist <- get_dba_file_list(rstr_data, measure, smooth, eddy)
  rstr_data@smooth <- smooth
  rstr_data@measure <- measure
  attrib_siz <- length(rstr_data@atlas_image)
  if ( !is.null(maskfile) ) {
    rstr_data@maskfile <- maskfile
    mask_image <- as.vector(RNifti::readNifti(maskfile))
    if ( length(mask_image) != attrib_siz) {
      stop(sprintf('Dimensions of atlas file %s and maskfile %s do not match', atlas_filename, maskfile), call. = FALSE)
    }
    rstr_data@mask_idx <- which(mask_image > 0)
  }
  else
    rstr_data@mask_idx = 1:attrib_siz
  
  first_subject_file <- as.vector(RNifti::readNifti(rstr_data@filelist[1]))
  # Check if the dimensions of first subject file and atlas match (the dimensions of atlas and mask are already checked above)
  if ( length(first_subject_file) != attrib_siz) {
    stop(sprintf('Dimensions of the atlas file %s and subject %s do not match', atlas_filename, rstr_data@filelist[1]), call. = FALSE)
  }
  
  rstr_data@data_array <- read_nii_images_for_all_subjects(rstr_data@filelist, attrib_siz, rstr_data@mask_idx)
  rstr_data@analysis_type <- "dba"
  rstr_data@data_type <- rs_data_types$nifti_image
  rstr_data@hemi <- "NA"
  return(rstr_data)
})

#' @rdname load_data
setMethod("load_data", signature = "RstrROIData", function(rstr_data, roiids = NULL, roimeas = NULL, exclude_col, roistat_fileprefix, roilabeldescfile) {
  
  rstr_data@analysis_type <- "roi"
  rstr_data@roiids <- roiids
  rstr_data@roimeas <- roimeas
  rstr_data@roilabeldescfile <- roilabeldescfile
  rstr_data@roistat_fileprefix <- roistat_fileprefix
  all_subjects <- rstr_load_roi_data(subjects_dir = rstr_data@subjdir,
                                     csv = rstr_data@csv,
                                     roiids = rstr_data@roiids,
                                     roimeas = rstr_data@roimeas,
                                     exclude_col=exclude_col,
                                     roistat_fileprefix=roistat_fileprefix, roilabeldescfile=roilabeldescfile)
  
  rstr_data@demographics <- as.data.frame(all_subjects[[1]])
  rstr_data@data_array <- matrix(nrow=nrow(rstr_data@demographics),ncol = length(rstr_data@roiids))
  for (col in 1:length(rstr_data@roiids)){
    current_col <- which(colnames(rstr_data@demographics) == paste0(rstr:::get_roi_tag(label_desc_df = rstr:::read_label_desc(roilabeldescfile),roiid=rstr_data@roiids[col])[[1]],"(",rstr_data@roiids[col],")"))
    rstr_data@data_array[,col] <- rstr_data@demographics[,current_col]
  }
  rstr_data@load_data_command <- sprintf("rstr_data <- load_rstr_data(type= 'roi',subjdir = '%s',csv= '%s',roiids= c( %s), roimeas= '%s', exclude_col='%s', roistat_fileprefix= '%s', roilabeldescfile='%s')",
                                         rstr_data@subjdir, rstr_data@csv, paste(rstr_data@roiids,collapse = ", "), rstr_data@roimeas, exclude_col, roistat_fileprefix, roilabeldescfile)
  
  
  return(rstr_data)
  
})


setMethod ("load_demographics", "RstrData", function(object) {
  object@demographics <- read_demographics(object@csv)
  return(object)
})



# signature(c(object = "RstrData", csv = "character"))

# setMethod("load_data", signature(object = "RstrROIData", subjects_dir = "character", csv = "character", roiids = "numeric",
#                                  roimeas = "character", outdir = "character"),
#           function(object, subjects_dir, csv, roiid, roimeas, outdir=NULL) {
#             print("hi")
#             }
#           )

#' Load data for statistical analysis.
#'
#' Loading data is usually the first step before running any statistical analysis.
#' Prior to using this function, MouseSuite and svreg should be run on all subjects.
#' If required, smoothing should be performed on cortical surface or volumetric image based measures.
#' A csv file containing subject demographic information should exist. The first column of this csv file
#' should have the subject identifiers. Subject identifiers can be alphanumeric
#' and should be exactly equal to the individual subject directory names.
#'
#' @param type character string denoting type of analysis. Should be sba, tbm, or roi.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param smooth numeric value denoting the smoothing level.
#' @param roiids numeric label identifiers for the regions of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. 
#' @param measure character specifying the brain imaging measure. If analyzing diffusion data, should be "FA".
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#' @param maskfile filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @param roistat_fileprefix file prefix for the ROI stats file (either .tissue.stats.txt or .roiwise.stats.txt)",
#' @param roilabeldescfile an xml file containing ROI label descriptions",
#' @examples
#' \dontrun{
#' my_sba_data <- load_rstr_data(type="sba", subjdir = "/path/to/my/subjectdirectory",
#' csv = "/path/to/my/demographics.csv", hemi = "left", smooth = 2.5)
#'
#' my_roi_data <- load_rstr_data(type="roi", subjdir = "/path/to/my/subjectdirectory",
#' csv="/path/to/my/demographics.csv", roiids=501, roimeas="gmthickness")
#' }
#'
#' @export
load_rstr_data <- function(type="sba", subjdir="", csv="", hemi="left",
                           smooth=0.0, roiids=0, roimeas="ROI Volume (mm^3)", measure="", atlas="", maskfile = "", eddy=TRUE, exclude_col = "",
                           roistat_fileprefix="roiwise.stats.txt", roilabeldescfile="") {
  
  atlas <- path.expand(atlas)
  maskfile <- path.expand(maskfile)
  subjdir <- path.expand(subjdir)
  
  valid_types <- c("sba", "tbm", "roi","dba","nca")
  if (! type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)
  
  switch(type,
         sba = { rstr_data <- load_sba_data_both_hemi(subjdir=subjdir, csv=csv, hemi=hemi, smooth = smooth, atlas=atlas, exclude_col=exclude_col) },
         tbm = { rstr_data <- load_tbm_data(subjdir=subjdir, csv=csv, smooth=smooth, atlas=atlas, maskfile=maskfile, exclude_col=exclude_col) },
         dba = { rstr_data <- load_dba_data(subjdir=subjdir, csv=csv, measure=measure, smooth=smooth, atlas=atlas, maskfile=maskfile, eddy=eddy, exclude_col=exclude_col) },
         roi = { rstr_data <- load_roi_data(subjdir, csv, roiids, roimeas, exclude_col=exclude_col, roistat_fileprefix=roistat_fileprefix, roilabeldescfile=roilabeldescfile) }
  )
  return(rstr_data)
}

#' Load data for statistical analysis from a filelist
#'
#' Loading data is usually the first step before running any statistical analysis.
#' Prior to using this function, MouseSuite and svreg should be run on all subjects.
#' If required, smoothing should be performed on cortical surface or volumetric image based measures.
#' Unlike \code{\link{load_rstr_data}}, this function loads data from a csv that contains a column for filelist
#' A csv file containing subject demographic information should exist. The first column of this csv file
#' should have the subject identifiers. Subject identifiers can be alphanumeric
#' and should be exactly equal to the individual subject directory names.
#' This csv file should also contain a column that contains a full path to the data file to be loaded.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param type character string denoting type of analysis. Should be sba, tbm, or roi.
#' @param file_col character string for the full file path.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param maskfile optional filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @examples
#' \dontrun{
#' my_data <- load_rstr_data_from_filelist(csv = "/path/to/my/demographics.csv",
#' type="sba", file_col = "COL_NAME", atlast = "/path/to/atlas",
#' maskfile = "/path/to/maskfile")
#' }
#'
#' @export
load_rstr_data_from_filelist <- function(csv="", subjdir="", hemi = "left", type="sba", file_col="", atlas="", maskfile = "") {
  
  #  rstr_sba_data <- new("RstrSBAData", subjdir=subjdir, csv=csv, exclude_col="")
  
  switch(type,
         tbm = { rstr_data <- load_tbm_data_from_filelist(subjdir = subjdir, csv = csv, file_col = file_col, atlas = atlas, maskfile = maskfile) },
         sba = { rstr_data <- load_sba_data_from_filelist(subjdir=subjdir, csv=csv, hemi = hemi, file_col = file_col, atlas=atlas) }
  )
  return(rstr_data)
}

#' Load cortical surface data for statistical analysis.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param smooth numeric value denoting the smoothing level.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
load_sba_data <- function(subjdir="", csv="", hemi="left", smooth=0.0, atlas="", exclude_col) {
  
  rstr_sba_data <- new("RstrSBAData", subjdir, csv, exclude_col)
  if (atlas == "") {
    mousesuite_atlas_id <- get_mousesuite_atlas_id_from_logfile(get_mousesuite_logfilename(subjdir, csv, exclude_col))
    sba_surf_atlas <- get_sba_atlas(mousesuite_atlas_id, hemi)
  }
  else
    sba_surf_atlas <- get_custom_sba_atlas_and_mask(atlas)
  
  rstr_sba_data <- load_data(rstr_sba_data, atlas_filename = sba_surf_atlas, hemi = hemi, smooth=smooth)
  rstr_sba_data@data_type <- rs_data_types$surface
  rstr_sba_data@hemi <- hemi
  return(rstr_sba_data)
}

#' Load cortical surface data for statistical analysis.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param file_col character string for the full file path.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @export
load_sba_data_from_filelist <- function(subjdir="", csv="", file_col="", atlas="", hemi="left") {
  
  rstr_data <- new("RstrSBAData", subjdir, csv, exclude_col="")
  
  rstr_data@atlas_filename <- atlas
  rstr_data@atlas_surface <- readdfs(atlas)
  sba_filelist <- rstr_data@demographics[, file_col]
  attrib_siz <- rstr_data@atlas_surface$hdr$nVertices
  rstr_data@data_array <- read_dfs_attributes_for_all_subjects(sba_filelist, attrib_siz)
  rstr_data@filelist <- sba_filelist
  rstr_data@analysis_type <- "sba"
  rstr_data@smooth <- 0.0
  rstr_data@data_type <- rs_data_types$surface
  rstr_data@hemi <- hemi
  return(rstr_data)
}

load_sba_data_both_hemi <- function(subjdir="", csv="", hemi="left", smooth=0.0, atlas="", exclude_col) {
  if (hemi == "both") {
    rstr_data_lh <- load_sba_data(subjdir=subjdir, csv=csv, hemi="left", smooth = smooth, atlas=atlas, exclude_col=exclude_col)
    rstr_data_rh <- load_sba_data(subjdir=subjdir, csv=csv, hemi="right", smooth = smooth, atlas=atlas, exclude_col=exclude_col)
    rstr_sba_data_both_hemi <- rstr_data_lh
    rstr_sba_data_both_hemi@data_array <- cbind(rstr_data_lh@data_array, rstr_data_rh@data_array)
    rstr_sba_data_both_hemi@nvertices_lh <- dim(rstr_data_lh@data_array)[2]
    rstr_sba_data_both_hemi@nvertices_rh <- dim(rstr_data_rh@data_array)[2]
    rstr_sba_data_both_hemi@atlas_surface_lh <- rstr_data_lh@atlas_surface
    rstr_sba_data_both_hemi@atlas_surface_rh <- rstr_data_rh@atlas_surface
    rstr_sba_data_both_hemi@hemi <- "both"
    rstr_sba_data_both_hemi@atlas_filename_lh <- rstr_data_lh@atlas_filename
    rstr_sba_data_both_hemi@atlas_filename_rh <- rstr_data_rh@atlas_filename
    return(rstr_sba_data_both_hemi)
  }
  else {
    rstr_data <- load_sba_data(subjdir=subjdir, csv=csv, hemi=hemi, smooth = smooth, atlas=atlas, exclude_col=exclude_col)
    return(rstr_data)
  }
}


#' Load tensor-based morphometry data for statistical analysis.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param smooth numeric value denoting the smoothing level.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param maskfile filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
load_tbm_data <- function(subjdir="", csv="", smooth=0.0, atlas="", maskfile="", exclude_col) {
  
  rstr_tbm_data <- new("RstrTBMData", subjdir, csv, exclude_col)
  tbm_atlas_and_mask = list()
  
  tbm_atlas_and_mask <- check_tbm_atlas_and_mask(subjdir, csv, atlas, maskfile, exclude_col)
  
  rstr_tbm_data <- load_data(rstr_tbm_data, atlas_filename = tbm_atlas_and_mask$nii_atlas, maskfile = tbm_atlas_and_mask$nii_atlas_mask, smooth=smooth)
  rstr_tbm_data@data_type <- rs_data_types$nifti_image
  return(rstr_tbm_data)
}

#' Load TBM data from file list for statistical analysis.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param file_col character string for the full file path.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' @param maskfile filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @export
#'
load_tbm_data_from_filelist <- function(subjdir="", csv="", file_col="", atlas="", maskfile="", exclude_col = "") {
  
  rstr_data <- new("RstrTBMData", subjdir, csv, exclude_col="")
  tbm_atlas_and_mask <- check_tbm_atlas_and_mask(subjdir, csv, atlas, maskfile, exclude_col)
  
  rstr_data@atlas_filename <- tbm_atlas_and_mask$nii_atlas
  rstr_data@atlas_image <- RNifti::readNifti(tbm_atlas_and_mask$nii_atlas)
  rstr_data@filelist <- rstr_data@demographics[, file_col]
  rstr_data@smooth <- 0.0
  rstr_data@measure <- ""
  
  attrib_siz <- length(rstr_data@atlas_image)
  if ( !is.null(maskfile) ) {
    rstr_data@maskfile <- maskfile
    mask_image <- as.vector(RNifti::readNifti(tbm_atlas_and_mask$nii_atlas))
    if ( length(mask_image) != attrib_siz) {
      stop(sprintf('Dimensions of atlas file %s and maskfile %s do not match', rstr_data@atlas_filename, maskfile), call. = FALSE)
    }
    rstr_data@mask_idx <- which(mask_image > 0)
  }
  else
    rstr_data@mask_idx = 1:attrib_siz
  first_subject_file <- as.vector(RNifti::readNifti(rstr_data@filelist[1]))
  # Check if the dimensions of first subject file and atlas match (the dimensions of atlas and mask are already checked above)
  if ( length(first_subject_file) != attrib_siz) {
    stop(sprintf('Dimensions of the atlas file %s and subject %s do not match. Check if you are using the correct atlas', rstr_data@atlas_filename, rstr_data@filelist[1]), call. = FALSE)
  }
  rstr_data@data_array <- read_nii_images_for_all_subjects(rstr_data@filelist, attrib_siz, rstr_data@mask_idx)
  rstr_data@analysis_type <- "tbm"
  rstr_data@data_type <- rs_data_types$nifti_image
  rstr_data@hemi <- "NA"
  
  return(rstr_data)
}

#' Load diffusion data for statistical analysis.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param measure character specifying the brain imaging measure. If analyzing diffusion data, should be "FA".
#' @param smooth numeric value denoting the smoothing level.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#' @param maskfile filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'

load_dba_data <- function(subjdir="", csv="", measure="", smooth=0.0, atlas="", eddy=TRUE, maskfile="", exclude_col) {
  
  rstr_dba_data <- new("RstrDBAData", subjdir, csv, exclude_col)
  dba_atlas_and_mask = list()
  
  dba_atlas_and_mask = list()
  if (maskfile == ""  && atlas == "") {
    mousesuite_atlas_id <- get_mousesuite_atlas_id_from_logfile(get_mousesuite_logfilename(subjdir, csv, exclude_col))
    dba_atlas_and_mask <- get_dba_atlas_and_mask(mousesuite_atlas_id)
  }
  else {
    if (maskfile != "" && atlas == "") {
      mousesuite_atlas_id <- get_mousesuite_atlas_id_from_logfile(get_mousesuite_logfilename(subjdir, csv, exclude_col))
      check_file_exists(maskfile, raise_error = TRUE)
      dba_atlas_and_mask$nii_atlas_mask <- maskfile
      dba_atlas_and_mask$nii_atlas <- get_dba_atlas(mousesuite_atlas_id)
    }
    else if (maskfile== "" && atlas != "") {
      mousesuite_atlas_id <- get_mousesuite_atlas_id_from_logfile(get_mousesuite_logfilename(subjdir, csv, exclude_col))
      check_file_exists(atlas, raise_error = TRUE)
      dba_atlas_and_mask$nii_atlas <- atlas
      dba_atlas_and_mask$nii_atlas_mask <- get_dba_mask(mousesuite_atlas_id)
    }
    else if (maskfile!= "" && atlas != "") {
      check_file_exists(maskfile, raise_error = TRUE)
      check_file_exists(atlas, raise_error = TRUE)
      dba_atlas_and_mask <- list("nii_atlas" = atlas, "nii_atlas_mask" = maskfile)
    }
  }
  
  
  rstr_dba_data <- load_data(rstr_dba_data, atlas_filename = dba_atlas_and_mask$nii_atlas, maskfile = dba_atlas_and_mask$nii_atlas_mask, measure=measure, smooth=smooth, eddy=eddy)
  rstr_dba_data@data_type <- rs_data_types$nifti_image
  return(rstr_dba_data)
}

#' Load ROI data for statistical analysis.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param roiids numeric label identifiers for the regions of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. 
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @param roistat_fileprefix file prefix for the ROI stats file (either .tissue.stats.txt or .roiwise.stats.txt)",
#' @param roilabeldescfile an xml file containing ROI label descriptions",
#'
load_roi_data <- function(subjdir="", csv="", roiids="", roimeas="", exclude_col, roistat_fileprefix, roilabeldescfile) {
  rstr_roi_data <- new("RstrROIData", subjdir, csv, exclude_col)
  rstr_roi_data <- load_data(rstr_roi_data, roiids = roiids, roimeas = roimeas, 
                             exclude_col=exclude_col, roistat_fileprefix=roistat_fileprefix, roilabeldescfile=roilabeldescfile)
  return(rstr_roi_data)
}

# #' Check that subject directory and demographics csv files exist
# #' @param object object of type \code{RstrData}
# #'
# check_files <- function(object){
#   if (!dir.exists(object@subjdir)) {
#     stop(sprintf("Subjects directory %s does not exist.\n", object@subjdir), call. = FALSE)
#   }
#
#   if (!file.exists(object@csv)) {
#     stop(sprintf("Demographics csv file %s does not exist.\n", object@csv), call. = FALSE)
#   }
# }

#' Package data for reproducible statistical analysis.
#'
#' Takes same parameters as \code{\link{load_rstr_data}} and copies the data to a new directory specified by outdir. You can repeatedly call this function to copy data of different types (tbm -- nii.gz, sba -- .dfs files etc.) to the same output directory.
#' Prior to using this function, MouseSuite and svreg should be run on all subjects.
#' If required, smoothing should be performed on cortical surface or volumetric image based measures.
#' A csv file containing subject demographic information should exist. The first column of this csv file
#' should have the subject identifiers. Subject identifiers can be alphanumeric
#' and should be exactly equal to the individual subject directory names.
#'
#' @param type character string denoting type of analysis. Should be sba, tbm, or roi.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param sbasmooth numeric value denoting the smoothing level for sba.
#' @param tbmsmooth numeric value denoting the smoothing level for tbm.
#' @param dbasmooth numeric value denoting the smoothing level for dba.
#' @param measure character specifying the brain imaging measure. If analyzing diffusion data, should be "FA".
#' @param atlas path name to the atlas
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#' @param outdir output directory that will contain the copied data
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
#' @export
package_data <- function(type="sba", subjdir=NULL, csv="", hemi="left",
                         sbasmooth=0.0, tbmsmooth=0.0,
                         dbasmooth=0.0, measure="FA", atlas="", eddy=TRUE, outdir=NULL, exclude_col="") {
  
  valid_types <- c("sba", "tbm", "roi","dba","nca", "all")
  if (! type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)
  
  if (is.null(subjdir))
    stop(sprintf("Subject directory is not specified."), call. = FALSE)
  
  if (is.null(outdir))
    stop(sprintf("Output directory is not specified."), call. = FALSE)
  
  if (outdir == "") {
    stop("Output directory is an empty string.", call.=FALSE)
  }
  if (dir.exists(outdir))
    stop("Output directory exists. Please specify a new output directory that does not exist.", call.=FALSE)
  
  switch(type,
         sba = { copy_sba_data(subjdir=subjdir, csv=csv, hemi=hemi, smooth = sbasmooth, outdir=outdir, exclude_col=exclude_col) },
         tbm = { copy_tbm_data(subjdir=subjdir, csv=csv, smooth=tbmsmooth, atlas=atlas, outdir=outdir, exclude_col=exclude_col) },
         dba = { copy_dba_data(subjdir=subjdir, csv=csv, measure=measure,
                               smooth=dbasmooth, atlas=atlas, eddy=eddy, outdir=outdir, exclude_col=exclude_col) },
         roi = { copy_roi_data(subjdir, csv, outdir=outdir, exclude_col=exclude_col) },
         all = {
           copy_sba_data(subjdir=subjdir, csv=csv, hemi="left", smooth = sbasmooth, outdir=outdir, exclude_col=exclude_col)
           copy_sba_data(subjdir=subjdir, csv=csv, hemi="right", smooth = sbasmooth, outdir=outdir, exclude_col=exclude_col)
           copy_tbm_data(subjdir=subjdir, csv=csv, smooth=tbmsmooth, atlas=atlas, outdir=outdir, exclude_col=exclude_col)
           copy_dba_data(subjdir=subjdir, csv=csv, measure=measure,
                         smooth=dbasmooth, atlas=atlas, eddy=eddy, outdir=outdir, exclude_col=exclude_col)
           copy_roi_data(subjdir, csv, outdir=outdir, exclude_col=exclude_col)
         }
  )
  # Copy the spreadsheet
  file.copy(csv, file.path(outdir, basename(csv)))
}

#' Package cortical surface data for reproducible statistical analysis.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param smooth numeric value denoting the smoothing level.
#' @param outdir output directory that will contain the copied data
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
copy_sba_data <- function(subjdir="", csv="", hemi="left", smooth=0.0, outdir, exclude_col="") {
  
  rstr_data <- new("RstrSBAData", subjdir, csv, exclude_col )
  logfilenames <- get_mousesuite_logfilename_for_all_subjects(subjdir, csv, exclude_col)
  mousesuite_atlas_id <- get_mousesuite_atlas_id_from_logfile(get_mousesuite_logfilename(subjdir, csv, exclude_col))
  sba_atlas_filename <- get_sba_atlas(mousesuite_atlas_id, hemi)
  sba_filelist <- get_sba_file_list(rstr_data, hemi, smooth)
  src_filelist <- c(sba_filelist, logfilenames, sba_atlas_filename)
  
  # Create subdirectories for subject IDs in outdir
  dir.create(file.path(outdir), showWarnings = FALSE)
  Vectorize(dir.create)(file.path(outdir, rstr_data@demographics$subjID), showWarnings = FALSE)
  dest_filelist <- file.path(outdir, rstr_data@demographics$subjID, basename(sba_filelist))
  dest_filelist <- c(dest_filelist, file.path(outdir, rstr_data@demographics$subjID, basename(logfilenames)),
                     file.path(outdir, basename(sba_atlas_filename)))
  # Copy files
  file_copy(src_filelist, dest_filelist, messg = "Copying sba data")
}

#' Package tensor-based morphometry data for reproducible statistical analysis.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param smooth numeric value denoting the smoothing level.
#' @param atlas path name to the atlas
#' @param outdir output directory that will contain the copied data
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
copy_tbm_data <- function(subjdir="", csv="", smooth=0.0, atlas, outdir, exclude_col="") {
  
  rstr_data <- new("RstrTBMData", subjdir, csv, exclude_col)
  logfilenames <- get_mousesuite_logfilename_for_all_subjects(subjdir, csv, exclude_col)
  mousesuite_atlas_id <- get_mousesuite_atlas_id_from_logfile(get_mousesuite_logfilename(subjdir, csv, exclude_col))
  tbm_atlas_mask_filename <- get_tbm_atlas_and_mask(mousesuite_atlas_id)
  tbm_filelist <- get_tbm_file_list(rstr_data, smooth)
  src_filelist <- c(tbm_filelist, logfilenames, tbm_atlas_mask_filename$nii_atlas, tbm_atlas_mask_filename$nii_atlas_mask)
  
  # Create subdirectories for subject IDs in outdir
  dir.create(file.path(outdir), showWarnings = FALSE)
  Vectorize(dir.create)(file.path(outdir, rstr_data@demographics$subjID), showWarnings = FALSE)
  dest_filelist <- file.path(outdir, rstr_data@demographics$subjID, basename(tbm_filelist))
  dest_filelist <- c(dest_filelist, file.path(outdir, rstr_data@demographics$subjID, basename(logfilenames)),
                     file.path(outdir, basename(tbm_atlas_mask_filename$nii_atlas)),
                     file.path(outdir, basename(tbm_atlas_mask_filename$nii_atlas_mask))
  )
  # Copy files
  message("Copying tbm data", appendLF = FALSE)
  file_copy(src_filelist, dest_filelist)
}

#' Package diffusion data for reproducible statistical analysis.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param measure character specifying the brain imaging measure. If analyzing diffusion data, should be "FA".
#' @param atlas path name to the atlas
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#' @param smooth numeric value denoting the smoothing level.
#' @param outdir output directory that will contain the copied data
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
copy_dba_data <- function(subjdir="", csv="", measure="FA", atlas="", eddy=TRUE, smooth=0.0, outdir, exclude_col="") {
  
  rstr_data <- new("RstrDBAData", subjdir, csv, exclude_col)
  logfilenames <- get_mousesuite_logfilename_for_all_subjects(subjdir, csv, exclude_col)
  mousesuite_atlas_id <- get_mousesuite_atlas_id_from_logfile(get_mousesuite_logfilename(subjdir, csv, exclude_col))
  dba_atlas_mask_filename <- get_dba_atlas_and_mask(mousesuite_atlas_id)
  
  message("Copying dba data", appendLF = FALSE)
  valid_diffusion_measures <- c('FA', 'MD', 'axial', 'radial', 'mADC', 'FRT_GFA')
  for (jj in valid_diffusion_measures) {
    dba_filelist <- get_dba_file_list(rstr_data, measure=jj, smooth = smooth, eddy = TRUE)
    src_filelist <- c(dba_filelist, logfilenames, dba_atlas_mask_filename$nii_atlas, dba_atlas_mask_filename$nii_atlas_mask)
    
    # Create subdirectories for subject IDs in outdir
    dir.create(file.path(outdir), showWarnings = FALSE)
    Vectorize(dir.create)(file.path(outdir, rstr_data@demographics$subjID), showWarnings = FALSE)
    dest_filelist <- file.path(outdir, rstr_data@demographics$subjID, basename(dba_filelist))
    dest_filelist <- c(dest_filelist, file.path(outdir, rstr_data@demographics$subjID, basename(logfilenames)),
                       file.path(outdir, basename(dba_atlas_mask_filename$nii_atlas)),
                       file.path(outdir, basename(dba_atlas_mask_filename$nii_atlas_mask))
    )
    # Copy files
    file_copy(src_filelist, dest_filelist)
  }
}

#' Package ROI data for reproducible statistical analysis.
#' @param subjdir subject directory containing MouseSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param outdir output directory that will contain the copied data
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @param roistat_fileprefix file prefix for the ROI stats file (either .tissue.stats.txt or .roiwise.stats.txt)",
#' @param roilabeldescfile an xml file containing ROI label descriptions",
#'
copy_roi_data <- function(subjdir="", csv="", outdir, exclude_col="", roistat_fileprefix="", roilabeldescfile="") {
  
  #TODO add extra arguments 
  rstr_data <- new("RstrROIData", subjdir, csv, exclude_col, roistat_fileprefix, roilabeldescfile)
  roiwise_file_list <- get_roi_file_list(rstr_data, roistat_fileprefix, roilabeldescfile)
  dest_filelist <- file.path(outdir, rstr_data@demographics$subjID, basename(roiwise_file_list))
  
  # Create subdirectories for subject IDs in outdir
  dir.create(file.path(outdir), showWarnings = FALSE)
  Vectorize(dir.create)(file.path(outdir, rstr_data@demographics$subjID), showWarnings = FALSE)
  
  # Copy files
  message("Copying ROI data", appendLF = FALSE)
  file_copy(roiwise_file_list, dest_filelist)
}

#' Copy files from the inputted source to the inputted destination
#' @param src_filelist list of source files
#' @param dest_filelist list of destination files
#' @param messg character string of displayed message
#' @param progress logical flag set TRUE to display progress bar in terminal
#'
file_copy <- function(src_filelist, dest_filelist, messg="Copying ", progress = TRUE) {
  
  if (length(src_filelist) != length(dest_filelist))
    stop(sprintf('The lengths of the source and the destination files do not match.'), call. = FALSE)
  
  message(messg, appendLF = FALSE)
  # If progress == TRUE, then display progressbar
  if (progress == TRUE) {
    pb <- txtProgressBar(max = length(src_filelist), style = 3)
    for (i in 1:length(src_filelist)) {
      file.copy(src_filelist[i], dest_filelist[i])
      setTxtProgressBar(pb, pb$getVal()+1)
    }
    close(pb)
  }
  else
    file.copy(src_filelist, dest_filelist)
}

