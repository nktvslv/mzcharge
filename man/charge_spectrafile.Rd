% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rfunctions.R
\name{charge_spectrafile}
\alias{charge_spectrafile}
\title{Centroid and assign charge to multiple profile spectra from mzML file}
\usage{
charge_spectrafile(
  filename,
  min_width,
  min_intensity,
  mz_tol,
  intensity_tol,
  max_charge,
  num_iso
)
}
\arguments{
\item{filename}{A path to mzML file.}

\item{min_width}{A size of sliding window for centroiding.}

\item{min_intensity}{Minimum intensity to consider for centroiding.}

\item{mz_tol}{Max mz deviation between isotopologues (ppm).}

\item{intensity_tol}{Max deviation between measured intensity and predicted
using averagine model.}

\item{max_charge}{Max charge to assign.}

\item{num_iso}{Max number of isotopologues to search for.}
}
\value{
A \code{data.frame}
}
\description{
\code{charge_spectrafile} reads profile spectra from mzML file and converts
to centroided spectra with assigned charge and position in isotopic envelope.
Uses functions for parallel processing if \href{https://cran.r-project.org/package=future.apply}{future.apply} is installed.
}
