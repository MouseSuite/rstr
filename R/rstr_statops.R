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

#' Perform analysis of variance (ANOVA) for brain imaging data.
#'
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and uses an
#' F-test to compare the full model including the \code{main_effect + covariates} with the reduced (null) model
#' that only includes the \code{covariates}.
#'
#' Slightly different from the standard R anova function, \code{rstr_anova} currently does not directly accept the
#' results from \code{lm_vec}. This could be accomodated in the future versions.

#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  rstr_data Object of type \code{\link{RstrData}}
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{rstr_p_adjust}} for valid values.
#' @seealso \code{\link{lm_vec}} for linear regression, \code{\link{rstr_ttest}} for independent sample and paired t-tests.
#'
#' @export

# param niter numeric variable for the number of iterations for permutations test. Will be ignored if mult_comp="fdr"
rstr_anova <- function(main_effect="", covariates="", rstr_data, mult_comp="fdr") { #, niter=5000
  
  if (class(rstr_data) == "RstrROIData") {
    return(rstr_roi_anova(main_effect = main_effect, covariates = covariates, rstr_data = rstr_data))
  }
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  rstr_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, rstr_data = rstr_data)
  rstr_lm_null <- lm_vec(main_effect = "", covariates = covariates, rstr_data = rstr_data)
  rstr_model <- anova_vec(rstr_lm_full, rstr_lm_null, rstr_data)
  
  switch(mult_comp,
         # perm={
         #   pvalue_and_nulldist <- maxTperm(main_effect = main_effect, covariates = covariates, rstr_data = rstr_data, niter)
         #   rstr_model@pvalues <- pvalue_and_nulldist[[1]]
         #   rstr_model@pvalues[is.nan(rstr_model@pvalues)] <- 1
         #   rstr_model@pvalues <- rstr_model@pvalues*rstr_model@tvalues_sign
         #   rstr_model@tvalues[abs(rstr_model@pvalues) >= 0.05] <- 0
         #   nulldist <- pvalue_and_nulldist[[2]]
         #   rstr_model@pvalues_adjusted <- perm_p_adjust(main_effect = main_effect, covariates = covariates, rstr_data, nulldist)
         #   rstr_model@tvalues_adjusted <- rstr_model@tvalues
         #   rstr_model@tvalues_adjusted[abs(rstr_model@pvalues_adjusted) >= 0.05] <- 0
         #   rstr_model@pvalues_adjusted <- rstr_model@pvalues_adjusted*rstr_model@tvalues_sign
         # },
         
         fdr={
           rstr_model@pvalues[is.nan(rstr_model@pvalues)] <- 1
           rstr_model@pvalues <- rstr_model@pvalues*rstr_model@tvalues_sign
           rstr_model@tvalues[abs(rstr_model@pvalues) >= 0.05] <- 0
           rstr_model@pvalues_adjusted <- rstr_p_adjust(rstr_model@pvalues, mult_comp)*rstr_model@tvalues_sign
           rstr_model@tvalues_adjusted <- rstr_model@tvalues
           rstr_model@tvalues_adjusted[abs(rstr_model@pvalues_adjusted) >= 0.05] <- 0
         }
  )
  
  message('Done.')
  return(rstr_model)
  
  
}

#' A vectorized version of analysis of variance (ANOVA).
#'
#' This function compares results of model fitting after \code{\link{lm_vec}}.
#' It accepts a full model and and a reduced model and compares them using an F-test.
#' For most scenarios, the user does not need to call this function directly. This function
#' will be called internally from \code{\link{rstr_anova}}
#'
#' @param rstr_lm_full An object of type \code{\link{RstrModel}} returned from \code{\link{lm_vec}}.
#' This is a full model including both the main effect and covariates.
#' @param rstr_lm_null An object of type \code{\link{RstrModel}} returned from \code{\link{lm_vec}}.
#' This is a null model including only the covariates.
#' @param  rstr_data Object of type \code{\link{RstrData}}
#'
#' @seealso \code{\link{rstr_anova}} for most commonly used function for ANOVA, \code{\link{lm_vec}} for vectorized linear regression, \code{\link{ttest_vec}} for
#' vectorized independent sample and paired t-tests.
#'
#' @export
anova_vec <- function(rstr_lm_full, rstr_lm_null, rstr_data) {
  
  N <- nrow(rstr_data@data_array)
  Fstat <- (rstr_lm_null@rss - rstr_lm_full@rss)/rstr_lm_full@rss * (N - rstr_lm_full@Npfull - 1)/(rstr_lm_full@Npfull - rstr_lm_null@Npnull)  # F statistic
  
  model_unique_idx <- which(rstr_lm_null@unique %in% rstr_lm_null@fullvars) + 1    # Add 1, because the first column in the design matrix is the intercept
  
  se_full_unique <- sqrt(diag(solve(t(rstr_lm_full@X_design_full) %*% rstr_lm_full@X_design_full)))[model_unique_idx] *
    sqrt(rstr_lm_full@rss / (N - rstr_lm_full@Npfull - 1))
  
  tvalues <- rstr_lm_full@beta_coeff[model_unique_idx, ]/(se_full_unique + .Machine$double.eps)
  pvalues <- 1 - pf(Fstat, rstr_lm_full@Npfull - rstr_lm_null@Npnull, N - rstr_lm_full@Npfull - 1)
  pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  tvalues_sign <- sign_tvalues(tvalues)
  
  rstr_model <- new("RstrModel", model_type="rstr_anova", main_effect = rstr_lm_full@main_effect, covariates = rstr_lm_full@covariates,
                    demographics = rstr_data@demographics, mspec_file="")
  rstr_model@pvalues <- pvalues
  rstr_model@tvalues <- tvalues
  rstr_model@tvalues_sign <- tvalues_sign
  rstr_model@se <- se_full_unique
  rstr_model@Fstat <- Fstat
  return(rstr_model)
}

#' Linear regression for brain imaging data.
#'
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and performs
#' a linear regression including \code{main_effect + covariates}.
#'
#' Slightly different from the standard R \code{lm} function, \code{rstr_lm} currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the \code{main_effect}
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  rstr_data Object of type \code{\link{RstrData}}
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{rstr_p_adjust}} for valid values.
#'
#' @export

# param  niter numeric variable for the number of iterations for permutations test. Will be ignored if mult_comp="fdr"
rstr_lm <- function(main_effect="", covariates="", rstr_data, mult_comp = "fdr") { #, niter=5000
  
  if (class(rstr_data) == "RstrROIData") {
    return(rstr_roi_anova(main_effect = main_effect, covariates = covariates, rstr_data = rstr_data))
  }
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  rstr_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, rstr_data = rstr_data)
  rstr_lm_null <- lm_vec(main_effect = "", covariates = covariates, rstr_data = rstr_data)
  rstr_model <- anova_vec(rstr_lm_full, rstr_lm_null, rstr_data)
  
  switch(mult_comp,
         # perm={
         #   pvalue_and_nulldist <- maxTperm(main_effect = main_effect, covariates = covariates, rstr_data = rstr_data, niter)
         #   rstr_model@pvalues <- pvalue_and_nulldist[[1]]
         #   rstr_model@pvalues[is.nan(rstr_model@pvalues)] <- 1
         #   rstr_model@pvalues <- rstr_model@pvalues*rstr_model@tvalues_sign
         #   #rstr_model@tvalues[abs(rstr_model@pvalues) >= 0.05] <- 0
         #   nulldist <- pvalue_and_nulldist[[2]]
         #   rstr_model@pvalues_adjusted <- perm_p_adjust(main_effect = main_effect, covariates = covariates, rstr_data, nulldist)
         #   rstr_model@tvalues_adjusted <- rstr_model@tvalues
         #   rstr_model@tvalues_adjusted[abs(rstr_model@pvalues_adjusted) >= 0.05] <- 0
         #   rstr_model@pvalues_adjusted <- rstr_model@pvalues_adjusted*rstr_model@tvalues_sign
         # },
         
         fdr={
           rstr_model@pvalues[is.nan(rstr_model@pvalues)] <- 1
           rstr_model@pvalues <- rstr_model@pvalues*rstr_model@tvalues_sign
           rstr_model@tvalues[abs(rstr_model@pvalues) >= 0.05] <- 0
           rstr_model@pvalues_adjusted <- rstr_p_adjust(rstr_model@pvalues, mult_comp)*rstr_model@tvalues_sign
           rstr_model@tvalues_adjusted <- rstr_model@tvalues
           rstr_model@tvalues_adjusted[abs(rstr_model@pvalues_adjusted) >= 0.05] <- 0
         }
  )
  
  message('Done.')
  return(rstr_model)
  
  # # Check the model type and call the appropriate method
  # Xtemp <- solve(t(rstr_model@X_design_full) %*% rstr_model@X_design_full) %*% t(rstr_model@X_design_full) # Pre Hat matrix
  # beta_full <- Xtemp %*% rstr_data@data_array  # beta coefficients
  # y_full <- rstr_model@X_design_full %*% beta_full  # Predicted response
  # RSS_full <- colSums((rstr_data@data_array - y_full)^2)
  #
  # Xtemp <- solve(t(rstr_model@X_design_null) %*% rstr_model@X_design_null) %*% t(rstr_model@X_design_null) # Pre Hat matrix
  # beta_null <- Xtemp %*% rstr_data@data_array  # beta coefficients
  # y_null <- rstr_model@X_design_null %*% beta_null  # Predicted response
  # RSS_null <- colSums((rstr_data@data_array - y_null)^2)
  #
  # N <- nrow(rstr_data@data_array)
  # Fstat <- (RSS_null - RSS_full)/RSS_full * (N - rstr_model@Npfull - 1)/(rstr_model@Npfull - rstr_model@Npnull)  # F statistic
  # model_unique_idx <- which(rstr_model@unique %in% rstr_model@fullvars) + 1    # Add 1, because the first column in the design matrix is the intercept
  #
  # se_full_unique <- sqrt(diag(solve(t(rstr_model@X_design_full) %*% rstr_model@X_design_full)))[model_unique_idx] *
  #   sqrt(RSS_full / (N - rstr_model@Npfull - 1))
  #
  # tvalue_sign <- (beta_full[model_unique_idx, ] + .Machine$double.eps)/(abs(beta_full[model_unique_idx, ]) + .Machine$double.eps)
  #
  # pvalues <- 1 - pf(Fstat, rstr_model@Npfull - rstr_model@Npnull, N - rstr_model@Npfull - 1)
  #
  # pvalues[is.nan(pvalues)] <- 1
  #
  # pvalues <- pvalues*tvalue_sign
  # tvalues <- beta_full[model_unique_idx, ]/(se_full_unique + .Machine$double.eps)
  # rstr_model@pvalues <- pvalues
  # rstr_model@tvalues <- tvalues
  # rstr_model@tvalues[abs(pvalues) >= 0.05] <- 0
  # rstr_model@pvalues_adjusted <- p.adjust(rstr_model@pvalues, 'BH')
  
}

#' Vectorized linear regression for brain imaging phenotypes.
#'
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and performs
#' a linear regression including \code{main_effect + covariates}.
#'
#' Slightly different from the standard R \code{lm} function, \code{lm_vec} currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the \code{main_effect}
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  rstr_data Object of type \code{\link{RstrData}}
#'
#' @export
lm_vec <- function(main_effect = "", covariates = "", rstr_data) {
  
  rstr_model <- new("RstrModel", model_type="rstr_lm", main_effect = main_effect, covariates = covariates,
                    demographics = rstr_data@demographics, mspec_file="")
  
  # lm_formula <- formula(sprintf('~ %s', paste(main_effect, '+', covariates)))
  lm_formula <- rstr_model@lm_formula
  
  # Fit model
  N <- dim(rstr_data@data_array)[1]
  Np <- length(unlist(strsplit(as.character(lm_formula)[2], '\\+')))
  
  X <- model.matrix(lm_formula, data = rstr_data@demographics)
  X_hat <- solve(t(X) %*% X) %*% t(X) # pre hat matrix
  beta_coeff <- X_hat %*% rstr_data@data_array  # beta coefficients
  Y <- X %*% beta_coeff  # predicted response
  rss <- colSums((rstr_data@data_array - Y)^2) # residual sum of squares
  
  if (main_effect == "")
    main_effect = "(Intercept)" # If main_effect is empty, return the parameters of the Intercept
  
  se <- sqrt(diag(solve(t(X) %*% X)))[[main_effect]] * sqrt(rss / (N-Np-1)) # standard error
  tvalues <- as.numeric(beta_coeff[main_effect, ]/(se + .Machine$double.eps)) # tvalues
  pvalues <- 2*pt(abs(tvalues), N-Np-1, lower.tail = FALSE) # pvalues
  residuals <- rstr_data@data_array - Y
  pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  rstr_model@pvalues <- pvalues
  rstr_model@tvalues <- tvalues
  rstr_model@beta_coeff <- beta_coeff
  rstr_model@rss <- rss
  rstr_model@residuals <- residuals
  return(rstr_model)
}


#' This function performs the ANOVA analysis for ROIs
#'
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  rstr_data Object of type \code{\link{RstrData}}
#'
#' @export

rstr_roi_anova <- function(main_effect="", covariates="", rstr_data=rstr_data) {
  # Check the model type and call the appropriate method
  rstr_model <- new("RstrModel", model_type="rstr_lm", main_effect = main_effect, covariates = covariates,
                    demographics = rstr_data@demographics, mspec_file="")
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  
  selected_col <- rep(NA, length(rstr_data@roiids))
  for (i in 1:length(rstr_data@roiids)){
    selected_col[i] <- paste0(as.character(get_roi_tag(read_label_desc(rstr_data@roilabeldescfile),rstr_data@roiids[i])), "(",rstr_data@roiids[i],")")
    rstr_data@demographics[,selected_col[i]]
  }
  
  cmd1 <- list()
  cmd2 <- list()
  cmd3 <- list()
  stats_commands <- list()
  
  for (i in 1:length(rstr_data@roiids)){
    cmd1[[i]] <- sprintf("lm_full_%s <- lm(%s, data = rstr_data@demographics)",as.character(rstr_data@roiids[i]),
                         paste('`',as.character(selected_col[i]),'`', ' ~ ', rstr_model@fullmodel, sep = ''))
    cmd2[[i]] <- sprintf("lm_null_%s <- lm(%s, data = rstr_data@demographics)",as.character(rstr_data@roiids[i]),
                         paste('`',as.character(selected_col[i]),'`', ' ~ ', rstr_model@nullmodel, sep = ''))
    cmd3[[i]] <- sprintf("pander::pander(anova(lm_full_%s, lm_null_%s))",as.character(rstr_data@roiids[i]),as.character(rstr_data@roiids[i]))
    stats_commands[[i]] <- c(cmd1[[i]], cmd2[[i]], cmd3[[i]])
  }
  
  for (i in 1:length(rstr_data@roiids)){
    for (cmd in stats_commands[[i]]) {
      eval(parse(text = cmd))
    }
  }
  rstr_model@stats_commands <- stats_commands
  
  rstr_model@load_data_command <- sprintf("rstr_model <- rstr_anova( main_effect = '%s', covariates = '%s', rstr_data = rstr_data) ",
                                          main_effect, covariates)
  
  return(rstr_model)
  
}


# rstr_roi_lm <- function(main_effect="", covariates="", rstr_data=rstr_data) {
#   # Check the model type and call the appropriate method
#   rstr_model <- new("RstrModel", model_type="rstr_lm", main_effect = main_effect, covariates = covariates,
#                    demographics = rstr_data@demographics, mspec_file="")
#   message('Running the statistical model. This may take a while...', appendLF = FALSE)
#
#   rstr_data@demographics[,as.character(rstr_data@roiids)]
#
#   cmd1 <- sprintf("lm_full_%s <- lm(%s, data = rstr_data@demographics)",as.character(rstr_data@roiids),
#                   paste('`',as.character(rstr_data@roiids),'`', ' ~ ', rstr_model@fullmodel, sep = ''))
#   cmd2 <- sprintf("lm_null_%s <- lm(%s, data = rstr_data@demographics)",as.character(rstr_data@roiids),
#                   paste('`',as.character(rstr_data@roiids),'`', ' ~ ', rstr_model@nullmodel, sep = ''))
#   cmd3 <- sprintf("pander::pander(anova(lm_full_%s, lm_null_%s))",as.character(rstr_data@roiids),as.character(rstr_data@roiids))
#
#   stats_commands <- c(cmd1, cmd2, cmd3)
#
#   for (cmd in stats_commands) {
#     eval(parse(text = cmd))
#   }
#   rstr_model@stats_commands <- stats_commands
#
#   return(rstr_model)
#
# }

#' Test for Correlation between a variable \code{corr_var} and a brain imaging phenotype.
#'
#' Test for correlation between a brain imaging phenotype (cortical thickness, determinant
#' of the jacobian matrix) and \code{corr_var} using the Pearson's product moment correlation
#' coefficient. The brain imaging phenotype is automatically selected from the type of \code{rstr_data}.
#' @param corr_var Character variable name. This should be present in the demographics csv file associated
#' with \code{rstr_data}.
#' @param  rstr_data Object of type \code{\link{RstrData}}.
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{rstr_p_adjust}} for valid values.
#' @details
#' \code{rstr_data} can be of the type "sba", "tbm", or "roi".
#'
#' @export
rstr_corr <- function(corr_var, rstr_data, mult_comp="fdr") {
  
  message('Running correlations...', appendLF = FALSE)
  rstr_model <- new("RstrModel", model_type="rstr_corr", corr_var = corr_var,
                    demographics = rstr_data@demographics, mspec_file="")
  
  corr_result <- corr_vec(rstr_data@data_array, rstr_data@demographics[[corr_var]])
  corr_coeff <- corr_result$corr_coeff
  rstr_model@tvalues <- corr_result$tvalues
  rstr_model@pvalues <- corr_result$pvalues
  tvalues_sign <- sign_tvalues(rstr_model@tvalues)
  rstr_model@tvalues_sign <- tvalues_sign
  
  rstr_model@pvalues <- sign(corr_coeff)*rstr_model@pvalues
  rstr_model@pvalues[abs(rstr_model@pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  rstr_model@pvalues[is.na(rstr_model@pvalues)] <- 1 # Set the p-values with the NA correlations to 1
  rstr_model@corr_values <- corr_coeff
  rstr_model@pvalues_adjusted <- rstr_p_adjust(rstr_model@pvalues, mult_comp)*rstr_model@tvalues_sign
  # rstr_model@pvalues_adjusted <- p.adjust(rstr_model@pvalues, 'BH')
  rstr_model@corr_values[abs(rstr_model@pvalues) >= 0.05] <- 0
  rstr_model@corr_values_masked_adjusted <- rstr_model@corr_values
  rstr_model@corr_values_masked_adjusted[abs(rstr_model@pvalues_adjusted) >= 0.05] <- 0
  rstr_model@tvalues_adjusted <- rstr_model@tvalues
  rstr_model@tvalues_adjusted[abs(rstr_model@pvalues_adjusted) >= 0.05] <- 0
  
  rstr_model@load_data_command <- sprintf("rstr_model <- rstr_corr(corr_var = '%s', rstr_data = rstr_data, mult_comp= '%s')",
                                          corr_var, mult_comp)
  
  message('Done.')
  return(rstr_model)
  
}

#' Vectorized correlation between a variable and brain imaging data.
#'
#' For most scenarios, the user does not need to call this function directly.
#' Instead call \code{\link{rstr_corr}} which calls this function internally.
#' @param X matrix of dimensions (\eqn{N x T}), where \eqn{N} = number of subjects and \eqn{T} = number of vertices/voxels.
#' @param Y vector of length \eqn{N}.
#' @export
corr_vec <- function(X, Y) {
  
  N <- length(Y)
  X_dev <- sweep(X, 2, colMeans(X))
  Y_dev <- Y - mean(Y)
  corr_coeff <- as.numeric((Y_dev %*% X_dev)/sqrt(colSums(X_dev^2)*sum(Y_dev^2)))
  tvalues <- corr_coeff * sqrt((N-2)/(1-corr_coeff^2 + .Machine$double.eps))
  pvalues <- 2*pt(abs(tvalues), N-2, lower.tail = FALSE)
  pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  return(list("tvalues"=tvalues, "pvalues"=pvalues, "corr_coeff"=corr_coeff))
}

#' T-test for for brain imaging phenotypes
#'
#' Perform independent sample and paired sample t-tests for differences between means of brain
#' imaging phenotypes for a categorical variable.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with \code{rstr_data}.
#' @param  rstr_data Object of type \code{\link{RstrData}}.
#' @param  paired logical; is TRUE if \code{group_var} contains matching (dependent) samples. The default value is \code{FALSE}.
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{rstr_p_adjust}} for valid values.
#' @details
#' The degrees of freedom are calculated using the Welch–Satterthwaite approximation by default.
#' \code{rstr_data} can be of the type "sba", "tbm", or "roi".
#'
#' @export
rstr_ttest <- function(group_var, rstr_data, paired = FALSE, mult_comp="fdr") {
  
  if (paired == FALSE)
    rstr_model <- new("RstrModel", model_type="unpairedttest", group_var = group_var,
                      demographics = rstr_data@demographics, mspec_file="")
  else
    rstr_model <- new("RstrModel", model_type="pairedttest", group_var = group_var,
                      demographics = rstr_data@demographics, mspec_file="")
  
  group1 <- levels(as.factor(rstr_data@demographics[[group_var]]))[1]
  group2 <- levels(as.factor(rstr_data@demographics[[group_var]]))[2]
  idx_group1 <- which(rstr_data@demographics[[group_var]] == group1)
  idx_group2 <- which(rstr_data@demographics[[group_var]] == group2)
  
  message(sprintf("The group variable is %s with two levels (%s, %s).", group_var, group1, group2), appendLF = TRUE)
  message(sprintf("Testing for the difference between the %s groups (%s - %s).", group_var, group1, group2), appendLF = TRUE)
  message('Running t-tests...', appendLF = FALSE)
  
  test_result <- ttest_vec(rstr_data@data_array[idx_group1,], rstr_data@data_array[idx_group2,], group_var, paired)
  
  rstr_model@pvalues <- test_result$pvalues
  rstr_model@tvalues <- test_result$tvalues
  rstr_model@tvalues_sign <- sign_tvalues(rstr_model@tvalues)
  
  rstr_model@pvalues[is.na(rstr_model@pvalues)] <- 1
  rstr_model@pvalues <- rstr_model@pvalues*rstr_model@tvalues_sign
  rstr_model@tvalues[abs(rstr_model@pvalues) >= 0.05] <- 0
  rstr_model@pvalues_adjusted <- rstr_p_adjust(rstr_model@pvalues, mult_comp)*rstr_model@tvalues_sign
  # rstr_model@pvalues_adjusted <- p.adjust(abs(rstr_model@pvalues), 'BH')
  rstr_model@tvalues_adjusted <- rstr_model@tvalues
  rstr_model@tvalues_adjusted[abs(rstr_model@pvalues_adjusted) >= 0.05] <- 0
  
  rstr_model@load_data_command <- sprintf("rstr_model <- rstr_ttest(group_var = '%s', rstr_data = rstr_data, paired = %s, mult_comp= '%s')",
                                          group_var, paired,mult_comp)
  
  message('Done.')
  return(rstr_model)
}

#' Vectorized t-test for for brain imaging phenotypes
#'
#' Perform independent sample and paired sample t-tests between two numerfor differences between means of brain
#' imaging phenotypes for a categorical variable.
#' For most scenarios, the user does not need to call this function directly. This function
#' will be called internally from \code{\link{rstr_ttest}}
#' @param X1 matrix of dimensions (\eqn{N1 x T}), where \eqn{N1} = number of subjects and \eqn{T} = number of vertices/voxels.
#' @param  X2 matrix of dimensions (\eqn{N2 x T}), where \eqn{N2} = number of subjects and \eqn{T} = number of vertices/voxels.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with \code{rstr_data}.
#' @param  paired logical; is TRUE if \code{group_var} contains matching (dependent) samples. The default value is \code{FALSE}.
#' @details
#' For an independent samples t-test \eqn{N1} not equal to \eqn{N2}.
#' For a dependent (paired) samples t-test, \eqn{N1 = N2}.
#' The degrees of freedom are calculated using the Welch–Satterthwaite approximation by default.
#' \code{rstr_data} can be of the type "sba", "tbm", or "roi".
#'
#' @export
ttest_vec <- function(X1, X2, group_var, paired=FALSE) {
  
  n1 <- dim(X1)[1]
  n2 <- dim(X2)[1]
  
  if (paired == TRUE) {
    D = X1 - X2
    D_mean <- colMeans(D)
    D_dev <- sweep(D, 2, D_mean)
    s1 <- sqrt(colSums(D_dev^2)/(n1-1))
    tvalues <- D_mean/(s1/sqrt(n1))
    pvalues <- 2*pt(abs(tvalues),n1-1, lower.tail = FALSE)
    pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  }
  else {
    X1_mean <- colMeans(X1)
    X2_mean <- colMeans(X2)
    
    X1_dev <- sweep(X1, 2, X1_mean)
    X2_dev <- sweep(X2, 2, X2_mean)
    
    SS1 <- colSums(X1_dev^2) # sum of squares of group 1
    SS2 <- colSums(X2_dev^2) # sum of squares of group 2
    
    s1_sq <- SS1/(n1-1)
    s2_sq <- SS2/(n2-1)
    
    se_diff <- sqrt(s1_sq/n1 + s2_sq/n2)
    tvalues <- (X1_mean - X2_mean)/(se_diff + .Machine$double.eps)
    # Calculate the degrees of freedom using the Welch–Satterthwaite approximation
    deg <- (s1_sq/n1 + s2_sq/n2)^2/(s1_sq^2/(n1^2*(n1-1)) + s2_sq^2/(n2^2*(n2-1)))
    pvalues <- 2*pt(abs(tvalues), deg, lower.tail = FALSE)
    pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  }
  return(list("tvalues"=tvalues, "pvalues"=pvalues))
}

#' Adjust p-values for multiple comparisons testing
#'
#' Perform multiple comparisons correction for mass univariate tests.
#' @param pvalues numeric vector of p values.
#' @param method character string specifying the method for correction. The default method is 'fdr'. Supported methods are
#' all types given in \code{\link{p.adjust.methods}}
#'
#' @export
rstr_p_adjust <- function(pvalues, method='fdr') {
  
  valid_methods <- c(p.adjust.methods) #, "perm"
  if ( !(method  %in% valid_methods) ) {
    warning(sprintf("%s is not a valid multiple comparisons method. Using no correction.", method), call. = FALSE)
    return(pvalues)
  }
  
  if (method %in% p.adjust.methods && method != "perm") {
    return(p.adjust(abs(pvalues), method))
  }
  
  # if (method == "perm") {
  #   tvalues_null <- rstr_perm()
  # }
  
  
}

# Calculate p-values for vanilla permutation test and determine null distribution for max-t permutation test.
#
# Perform permutation tests using Freedman-Lane method.
# @param main_effect Character string containing an independent variable whose effect you want to measure.
# It could be disease status, age, gender etc. This should strictly be a single variable. This can be
# either a categorical or a continuous variable.
# @param covariates Character string containing a set of other predictors (variables) in the model. If more than
# one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
# @param rstr_data Object of type \code{\link{RstrData}}
# @param num_of_perm Number of iterations/shuffles for permutation test.
# @details
# The permutation test handles the exchangeability assumption with Freedman-Lane Method. This function also
# utilizes multiprocessing to reduce computing time.
# \code{rstr_data} can be of the type "sba", "tbm", or "roi".
#
# @export
# maxTperm <- function(main_effect = "", covariates = "", rstr_data, num_of_perm){
#
#   N <- dim(rstr_data@data_array)[1]
#
#   rstr_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, rstr_data = rstr_data)
#   rstr_lm_null <- lm_vec(main_effect = "", covariates = covariates, rstr_data = rstr_data)
#
#   # (1) compute full model t-statistic
#   T0 <- rstr_lm_full@tvalues
#   maxT0 <- T0[ which.max( abs(T0) ) ]
#
#   # (2) compute estimated gamma_hat and estimated residuals from reduced model
#   # gamma_hat <- rstr_lm_null@beta_coeff
#   # residuals_null <- rstr_lm_null@residuals
#
#   # (3) compute a set of permuted data Y
#   rstr_data_cp <- rstr_data
#
#   t_bin_int <- rep(0, dim(rstr_data@data_array)[2])
#   t_max_per_perm <- rep(0, num_of_perm-1)
#   # rewrote into for loop, single core
#   for (j in 1:(num_of_perm-1)){
#     set.seed(j)
#     pmatrix <- as(sample(N), "pMatrix")
#     Y_j <- (pmatrix %*% rstr_lm_null@residuals) + (rstr_lm_null@X_design_null %*% rstr_lm_null@beta_coeff)
#
#     # (4) regress permuted data Y_j against the full model
#     rstr_data_cp@data_array <- Y_j
#     rstr_lm_full_perm <- lm_vec(main_effect = main_effect, covariates = covariates, rstr_data = rstr_data_cp )
#
#     ## binarize vector after comparing permuted T and observed T, where T0 is the vector containing
#     # all t values from the observed model
#     idx <- which( abs(T0) <= abs(rstr_lm_full_perm@tvalues) )
#
#     t_bin_int[idx] <- t_bin_int[idx] + 1
#
#     # t_bin <- as.bit(rep(FALSE, dim(rstr_data_cp@data_array)[2]))
#     # t_bin[idx] <- TRUE
#
#     t_max_per_perm[j] <- rstr_lm_full_perm@tvalues[ which.max( abs(rstr_lm_full_perm@tvalues) ) ] # max tvalue from perm
#
#   }
#
#   pvalues <- rep(0, dim(rstr_data@data_array)[2])
#   for (n in 1:dim(rstr_data@data_array)[2]){
#     count <- t_bin_int[n]
#     pvalues[n] <- as.double((count+1)/num_of_perm)
#   }
#
#   return(list(pvalues, t_max_per_perm))
#
# }


# Adjust p-values for multiple comparisons testing
#
# Perform multiple comparisons correction for mass univariate tests using the max-t method.
# @param main_effect Character string containing an independent variable whose effect you want to measure.
# It could be disease status, age, gender etc. This should strictly be a single variable. This can be
# either a categorical or a continuous variable.
# @param covariates Character string containing a set of other predictors (variables) in the model. If more than
# one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
# @param rstr_data Object of type \code{\link{RstrData}}
# @param tvalues_null Null distribution output from \code{\link{maxTperm}}. Statistics drawn from each
# shuffle are the maximum t-statistic within ROI.
#
# @export
# perm_p_adjust <- function(main_effect = "", covariates = "", rstr_data, tvalues_null){
#   rstr_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, rstr_data = rstr_data)
#
#   tvalues <- rstr_lm_full@tvalues
#   pvalues_adj <- rep(0, length(tvalues))
#
#   for (i in 1:length(tvalues)){
#     p <- (sum(abs(tvalues_null) >= abs(tvalues[i]))+1) / (length(tvalues_null)+1)
#     pvalues_adj[i] <- p
#   }
#   return(pvalues_adj)
# }

#' Linear mixed-effects model for brain imaging data.
#'
#' Linear regression for brain imaging data.
#'
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and performs
#' a linear regression including \code{main_effect + covariates}.
#'
#' Slightly different from the standard R \code{lmer} function, \code{rstr_lmer} currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the \code{main_effect}
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with \code{rstr_data}.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  rstr_data Object of type \code{\link{RstrData}}
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{rstr_p_adjust}} for valid values.
#'
#' @export
rstr_lmer <- function(group_var, main_effect="", covariates="", rstr_data, mult_comp = "fdr") {
  
  if (class(rstr_data) == "RstrROIData") {
    return(rstr_roi_lmer_anova(group_var, main_effect = main_effect, covariates = covariates, rstr_data = rstr_data))
  }
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  rstr_model <- lmer_vec(group_var, main_effect = main_effect, covariates = covariates, rstr_data = rstr_data)
  # switch(mult_comp,
  #        perm={
  #          cl <- parallel::makeCluster(parallel::detectCores())
  #          registerDoParallel(cl)
  #          options(warn=-1)
  #          pvalue_and_nulldist <- maxTperm(main_effect = main_effect, covariates = covariates, rstr_data = rstr_data, niter)
  #          rstr_model@pvalues <- pvalue_and_nulldist[[1]]
  #          rstr_model@pvalues[is.nan(rstr_model@pvalues)] <- 1
  #          rstr_model@pvalues <- rstr_model@pvalues*rstr_model@tvalues_sign
  #         #rstr_model@tvalues[abs(rstr_model@pvalues) >= 0.05] <- 0
  #          nulldist <- pvalue_and_nulldist[[2]]
  #          rstr_model@pvalues_adjusted <- perm_p_adjust(main_effect = main_effect, covariates = covariates, rstr_data, nulldist)
  #          rstr_model@tvalues_adjusted <- rstr_model@tvalues
  #          rstr_model@tvalues_adjusted[abs(rstr_model@pvalues_adjusted) >= 0.05] <- 0
  #          stopCluster(cl)
  #          options(warn=0)
  #        },
  #        fdr={
  #          rstr_model@pvalues[is.nan(rstr_model@pvalues)] <- 1
  #          rstr_model@pvalues <- rstr_model@pvalues*rstr_model@tvalues_sign
  #          #rstr_model@tvalues[abs(rstr_model@pvalues) >= 0.05] <- 0
  #          rstr_model@pvalues_adjusted <- rstr_p_adjust(rstr_model@pvalues, mult_comp)
  #          rstr_model@tvalues_adjusted <- rstr_model@tvalues
  #          rstr_model@tvalues_adjusted[abs(rstr_model@pvalues_adjusted) >= 0.05] <- 0
  #        }
  # )
  
  message('Done.')
  return(rstr_model)
}

#' Vectorized linear mixed-effects regression for brain imaging phenotypes.
#'
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and performs
#' a linear regression including \code{main_effect + covariates}.
#'
#' Slightly different from the standard R \code{lmer} function, \code{lmer_vec} currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the \code{main_effect}
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with \code{rstr_data}.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  rstr_data Object of type \code{\link{RstrData}}
#'
#' @export
lmer_vec <- function(group_var, main_effect = "", covariates = "", rstr_data) {
  
  rstr_model <- new("RstrModel", model_type="rstr_lmer", main_effect = main_effect, covariates = covariates,
                    group_var = group_var, demographics = rstr_data@demographics, mspec_file="")
  
  # Mass Univariate Linear Mixed Effects Analysis (still in progress)
  #lm_formula <- formula(sprintf('~ %s', paste(main_effect, '+', covariates)))
  #lmer_formula <- rstr_model@lm_formula
  
  # Fit model
  N <- dim(rstr_data@data_array)[1]
  # Np <- length(unlist(strsplit(as.character(covariates), '\\+'))) + 1
  #
  # X <- model.matrix(lm_formula, data = rstr_data@demographics)
  # X_hat <- solve(t(X) %*% X) %*% t(X) # pre hat matrix
  # beta_coeff <- X_hat %*% rstr_data@data_array  # beta coefficients
  # Y <- X %*% beta_coeff  # predicted response
  # rss <- colSums((rstr_data@data_array - Y)^2) # residual sum of squares
  #
  # if (main_effect == "")
  #   main_effect = "(Intercept)" # If main_effect is empty, return the parameters of the Intercept
  #
  # se <- sqrt(diag(solve(t(X) %*% X)))[[main_effect]] * sqrt(rss / (N-Np-1)) # standard error
  # tvalues <- as.numeric(beta_coeff[main_effect, ]/(se + .Machine$double.eps)) # tvalue
  # pvalues <- 2*pt(abs(tvalues), N-Np-1, lower.tail = FALSE) # pvalu
  # residuals <- rstr_data@data_array - Y
  # pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  # rstr_model@pvalues <- pvalues
  # rstr_model@tvalues <- tvalues
  # rstr_model@beta_coeff <- beta_coeff
  # rstr_model@rss <- rss
  # rstr_model@residuals <- residuals
  
  return(rstr_model)
}

#' This function performs the ANOVA Linear Mixed Effects analysis for ROIs
#'
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with \code{rstr_data}.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  rstr_data Object of type \code{\link{RstrData}}
#'
#' @export

rstr_roi_lmer_anova <- function(group_var, main_effect="", covariates="", rstr_data){
  
  rstr_model <- new("RstrModel", model_type="rstr_lmer", main_effect = main_effect, covariates = covariates,
                    group_var = group_var, demographics = rstr_data@demographics, mspec_file="")
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  
  # For univariate lmer
  selected_col <- rep(NA, length(rstr_data@roiids))
  for (i in 1:length(rstr_data@roiids)){
    selected_col[i] <- paste0(as.character(get_roi_tag(read_label_desc(rstr_data@roilabeldescfile),rstr_data@roiids[i])), "(",rstr_data@roiids[i],")")
  }
  
  cmd1 <- list()
  cmd2 <- list()
  cmd3 <- list()
  stats_commands <- list()
  
  for (i in 1:length(rstr_data@roiids)){
    rstr_lmer_full_formula <- paste0("`",as.character(get_roi_tag(read_label_desc(rstr_data@roilabeldescfile),rstr_data@roiids[i])), "(",rstr_data@roiids[i],")`"," ~ ",rstr_model@main_effect," + (",rstr_model@main_effect,"|",rstr_model@group_var,") + ",rstr_model@covariates)
    rstr_lmer_null_formula <- paste0("`",as.character(get_roi_tag(read_label_desc(rstr_data@roilabeldescfile),rstr_data@roiids[i])), "(",rstr_data@roiids[i],")`"," ~ (1|",rstr_model@group_var,") + ",rstr_model@covariates)
    
    cmd1[[i]] <- sprintf("lmer_full_%s <- lme4::lmer(%s, data = rstr_data@demographics,control=lme4::lmerControl(check.nobs.vs.nRE='ignore'))",
                         as.character(rstr_data@roiids[i]),rstr_lmer_full_formula)
    cmd2[[i]] <- sprintf("lmer_null_%s <- lme4::lmer(%s, data = rstr_data@demographics,control=lme4::lmerControl(check.nobs.vs.nRE='ignore'))",
                         as.character(rstr_data@roiids[i]),rstr_lmer_null_formula)
    cmd3[[i]] <- sprintf("pander::pander(anova(lmer_full_%s, lmer_null_%s))",as.character(rstr_data@roiids[i]),as.character(rstr_data@roiids[i]))
    stats_commands[[i]] <- c(cmd1[[i]], cmd2[[i]], cmd3[[i]])
  }
  
  eval(parse(text = stats_commands))
  rstr_model@stats_commands <- stats_commands
  rstr_model@load_data_command <- sprintf("rstr_model <- rstr_lmer(group_var = '%s', main_effect = '%s', covariates = '%s', rstr_data = rstr_data) ",
                                          group_var, main_effect, covariates)
  
  message('Done.')
  return (rstr_model)
}
