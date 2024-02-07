#' Nadaraya-Watson type kernel estimation of the conditional mean and covariance
#'
#' Contains functions to estimate the conditional mean and covariance using a 
#' Nadaraya-Watson kernel based estimator. This package contains the 
#' Nadaraya-Watson kernel-based estimator for the conditional mean and covariance 
#' presented in Yin et al. (2010) and Neumann et al. (2024) which can be used to 
#' condition the mean or covariance between system outputs on confounders, e.g. 
#' temperature.
#' 
#'  Neumann et al. (2024).
#'  Confounder-adjusted Covariances of System Outputs and Applications to Structural Health Monitoring.
#'  \emph{PREPRINT}, \strong{XX}, pp. 1--26.
#'
#'  Yin et al. (2010).
#'  Nonparametric Covariance Model. \emph{Statistica Sinica} \strong{20}, 469-–79.
#' 
#' @docType package
#' @importFrom Rcpp evalCpp
#' @useDynLib covest, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#' @name covest-package
NULL
