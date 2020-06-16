#' \code{mzcharge} package
#'
#' Deconvolution of high-resolution MS1 spectra for label-free peptide quantitation.
#'
#' @docType package
#' @name mzcharge
#'
#' @useDynLib mzcharge
#' @importFrom Rcpp sourceCpp
NULL

utils::globalVariables(c(".", "mz", "mz_tol", "mzdif", "mzclu", "query_mz",
                         "charge", "query_charge", "scdif", "scclu", "scan_tol",
                         "retention_time", "acquisition_num", "scan_order",
                         "intensity", "iso_position", "min_iso", "iso_sum",
                         "iso_count", "mass", "msLevel"))
