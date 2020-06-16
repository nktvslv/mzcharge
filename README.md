# mzcharge
## Deconvolution of high-resolution MS1 spectra
`mzcharge` is an R package for centroiding and charge assignment for peaks in
high-resolution mass spectra of peptides. It is intended to be used for
label-free quantitation of peptides analyzed by LC-MS.

## Installation

### Dependencies
```R
# CRAN packages
install.packages(c("Rcpp","RcppArmadillo","data.table","future","future.apply"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("mzR")
```
When on Windows, install Rtools from [CRAN](https://cran.r-project.org) to 
compile C++ code

### Install package with `devtools`
```R
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")
devtools::install_github("https://github.com/nktvslv/mzcharge.git")
```

## Package functionality
`mzcharge` has five functions:

`charge_singlespec` is the major function that takes a single spectrum (numeric
vector of mz values followed by intensity values, an output of `spectra`
function from [mzR](https://github.com/sneumann/mzR/) package) and returns a
data.frame where each row represents a centroid with corresponding mz, 
intensity, charge, isotopic position, index of isotopic pattern and count of 
peaks in pattern. This function can be used interactively for tweaking processing
parameters (peak width, mz tolerance, averagine correlation).

`charge_spectrafile` runs `charge_singlespec` on every MS1 spectrum from
raw file in mzML format. It uses [mzR](https://github.com/sneumann/mzR/) to read 
mzML and [future.apply](https://cran.r-project.org/package=future.apply) for 
parallel processing. See [ProteoWizard](http://proteowizard.sourceforge.net)
for mass spectrometry raw data format conversion.

`charge_corr` corrects charges for monoisotopic peaks by comparing centroids
from neighboring spectra. Takes data.frame output from `charge_spectrafile`
as its input. This function uses 
[data.table](https://cran.r-project.org/package=data.table) package to spead
up processing.

`xic_monomass` and `xic_monomz` extract ion current for either defined 
monoisotopic mass (uncharged) or a pair of monoisotopic mz and charge values. 
Both functions implemented as data.table subsetting. The data.table output can
be directed to plotting functions.
