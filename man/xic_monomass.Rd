% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rfunctions.R
\name{xic_monomass}
\alias{xic_monomass}
\title{Extract ion current for uncharged monoisotopic mass}
\usage{
xic_monomass(data, query_mass, mass_tol)
}
\arguments{
\item{data}{a \code{data.table} output from \code{\link[mzcharge]{charge_corr}}.}

\item{query_mass}{query mass.}

\item{mass_tol}{max allowed mass deviation from \code{query_mass} in ppm.}
}
\value{
a \code{data.table}.
}
\description{
Extract ion current for uncharged monoisotopic mass
}
