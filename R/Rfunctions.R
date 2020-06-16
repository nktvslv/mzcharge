#' Centroid and assign charge to multiple profile spectra from mzML file
#'
#' @description \code{charge_spectrafile} reads profile spectra from mzML file and converts
#' to centroided spectra with assigned charge and position in isotopic envelope.
#' Uses functions for parallel processing if \href{https://cran.r-project.org/package=future.apply}{future.apply} is installed.
#'
#' @param filename A path to mzML file.
#' @param min_width A size of sliding window for centroiding.
#' @param min_intensity Minimum intensity to consider for centroiding.
#' @param mz_tol Max mz deviation between isotopologues (ppm).
#' @param intensity_tol Max deviation between measured intensity and predicted
#' using averagine model.
#' @param max_charge Max charge to assign.
#' @param num_iso Max number of isotopologues to search for.
#'
#' @return A \code{data.frame}
#' @export
#' @import data.table
#' @import future
#' @import future.apply
charge_spectrafile <- function(filename, min_width, min_intensity,
                               mz_tol, intensity_tol, max_charge, num_iso) {

    # use mzR to read file and extract MS1 spectra
    mzml <- mzR::openMSfile(filename)
    ms1headers <- subset(mzR::header(mzml), msLevel == 1)

    #check if spectra are already centroided
    if (all(ms1headers[["centroided"]])) {
        stop("Spectra in", filename, "appears to be centroided.")
    }
    ms1spectra <- mzR::spectra(mzml, ms1headers[["acquisitionNum"]])
    mzR::close(mzml)

    plan(multisession)

    out <- future_lapply(X = ms1spectra, FUN = charge_singlespec,
                         min_width, min_intensity, mz_tol, intensity_tol, max_charge, num_iso)

    out <- mapply(FUN = function(a, b, c, d) cbind(a, retention_time = b / 60, acquisition_num = c, scan_order = d),
                  out, ms1headers[["retentionTime"]], ms1headers[["acquisitionNum"]], 1:length(out),
                  SIMPLIFY = FALSE)

    out <- rbindlist(out)

    return(out)
}

#' Correct misassigned and unassigned charges for monoisotopic peaks
#'
#' @param data a \code{data.table} output of \code{\link[mzcharge]{charge_spectrafile}} function.
#' @param mz_tol max allowed deviation in m/z between centroids in neighboring scans (ppm).
#' @param scan_tol max allowed gap between scans to consider them nighbors.
#' @param min_iso min number of isotopologues used for initial charge assignment.
#'
#' @return a \code{data.table} containing monoisotopic peaks with corrected charge.
#' @export
#'
#' @import data.table
charge_corr <- function(data, mz_tol, scan_tol, min_iso) {
    # fast proximity clustering by mz and rt
    setorder(data, mz)
    data[,mzdif := c(0, diff(mz)) / mz * 1e6][,mzclu := ifelse(mzdif < sqrt(mz_tol), 0L, 1L)][,mzclu := cumsum(mzclu)]
    setorder(data, mzclu, scan_order)
    data[,mzdif := c(0, diff(mz)) / mz * 1e6][,mzclu := ifelse(mzdif < mz_tol, 0L, 1L)][,mzclu := cumsum(mzclu)]
    data[,scdif := c(0, diff(scan_order)-1)][,scclu := ifelse(scdif < scan_tol, 0L, 1L)][,scclu := cumsum(scclu)]

    # calculate the sum of isotopic peaks per charge per mz per rt
    data[charge > 0 & iso_position == 0 & iso_count >= min_iso,
         iso_sum := sum(iso_count),
         by = .(charge, mzclu, scclu)]

    # assign charge which has max sum of isotopic peaks
    data[iso_position == 0,
         charge := charge[which.max(iso_sum)],
         by = .(mzclu, scclu)]

    setorder(data, retention_time, mz)

    return(data[iso_position == 0 & charge > 0,
                .(mz, charge, retention_time, acquisition_num, scan_order, intensity)])
}


#' Extract ion current for monoisotopic mz/charge pair
#'
#' @param data a \code{data.table} output of \code{\link[mzcharge]{charge_corr}} function.
#' @param query_mz query mz.
#' @param query_charge query charge.
#' @param mz_tol max allowed deviation from \code{query_mz} in ppm.
#'
#' @return a \code{data.table}.
#' @export
#'
#' @import data.table
xic_monomz <- function(data, query_mz, query_charge, mz_tol) {
    data[charge == query_charge &
             abs(mz - query_mz) / query_mz * 1e6 < mz_tol]
}


#' Extract ion current for uncharged monoisotopic mass
#'
#' @param data a \code{data.table} output from \code{\link[mzcharge]{charge_corr}}.
#' @param query_mass query mass.
#' @param mass_tol max allowed mass deviation from \code{query_mass} in ppm.
#'
#' @return a \code{data.table}.
#' @export
#'
#' @import data.table
xic_monomass <- function(data, query_mass, mass_tol) {
    # convert mz * charge to mass and xic
    PROTON_MASS <- 1.007276
    data[,mass := mz * charge - PROTON_MASS * charge][abs(mass - query_mass)/ query_mass * 1e6 < mass_tol]
}
