#include <RcppArmadillo.h>
#include "constants.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(c++11)]]

double ppm_diff(double x, double y) {
    return std::abs(x - y) / y * 1e6;
}

double frac_diff(double x, double y) {
    return std::abs(x - y) / y;
}

// [[Rcpp::export]]
List find_centroids(const NumericVector& spectrum, int min_width,
                    double min_intensity) {
    // check parameters
    if (spectrum.size() % 2 != 0) {
      stop("spectrum appears to have different number of mz and intensity values");
    }

    if ((min_width < 3) || (min_width % 2 != 1)) {
      stop("min_width should be odd number above 3");
    }

    if (min_intensity < 0) {
      stop("min_intensity should be zero or positive number");
    }

    // separate spectrum into mz (x) and intensity (y) vectors
    NumericVector x = NumericVector(spectrum.cbegin(), spectrum.cbegin() + spectrum.size() / 2);
    NumericVector y = NumericVector(spectrum.cbegin() + spectrum.size() / 2, spectrum.cend());

    // scan through intensity vector (y) and find local maxima accross min_width points
    // tested with Orbitrap-recorded spectrum spectra
    NumericVector xc, yc;
    int span = (int) min_width / 2 + 1;
    for (int i = span; i < y.size() - span; ++i) {
        // goto next i if i-th point is below intensity threshold
        if (y[i] < min_intensity) {
            continue;
        }

        bool is_centroid = true;
        // make sure span-long ascend and descend at i-th peak are continuous
        for (int j = 0; j < span; ++j) {
            if ((y[i-j] <= y[i-j-1]) || (y[i+j] <= y[i+j+1])) {
                is_centroid = false;
                break;
            }
        }

        if (is_centroid) {
            // weighted mean for mz, max for intensity
            xc.push_back((x[i]*y[i]+x[i-1]*y[i-1]+x[i+1]*y[i+1])/(y[i]+y[i-1]+y[i+1]));
            yc.push_back(y[i]);
            i += min_width;
        }
    }
    return List::create(_["mz"] = xc, _["intensity"] = yc);
}

NumericVector averagine_ratios(double mz, int charge, int num_iso) {

    // peptide mass
    double mass = mz * charge - PROTON_MASS * charge;

    // elemental composition based on averagine amino acid model
    double pep_len = mass / AVERAGINE_MASS;
    double nC = pep_len * N_CARBON;
    double nH = pep_len * N_HYDROGEN;
    double nN = pep_len * N_NITRO;
    double nO = pep_len * N_OXIGEN;
    double nS = pep_len * N_SULFUR;

    // initialize vectors with isotopes distribution for each element and calculate FFT
    // implementation from Sadygov RG. 2018. J. Proteome Res. 17: 751–758.
    arma::vec c12(num_iso, arma::fill::zeros);
    c12[0] = P_C12; c12[1] = P_C13;
    arma::cx_vec c12ri = arma::pow(arma::fft(c12), nC);

    arma::vec h1(num_iso, arma::fill::zeros);
    h1[0] = P_H1; h1[1] = P_H2;
    arma::cx_vec h1ri = arma::pow(arma::fft(h1), nH);

    arma::vec n14(num_iso, arma::fill::zeros);
    n14[0] = P_N14; n14[1] = P_N15;
    arma::cx_vec n14ri = arma::pow(arma::fft(n14), nN);

    arma::vec o16(num_iso, arma::fill::zeros);
    o16[0] = P_O16; o16[1] = P_O17; o16[2] = P_O18;
    arma::cx_vec o16ri = arma::pow(arma::fft(o16), nO);

    arma::vec s32(num_iso, arma::fill::zeros);
    s32[0] = P_S32; s32[1] = P_S33; s32[2] = P_S34; s32[4] = P_S36;
    arma::cx_vec s32ri = arma::pow(arma::fft(s32), nS);

    // combine FFTs of individual elements
    using arma::operator%;
    arma::cx_vec prod = c12ri % h1ri;
                 prod = prod % n14ri;
                 prod = prod % o16ri;
                 prod = prod % s32ri;

    // inverse FFT
    arma::cx_vec inv_prod = arma::ifft(prod);

    // get real part of inverse FFT
    arma::vec results = arma::real(inv_prod);

    // normalize ratios
    results = results / results[0];

    return wrap(results);
}

IntegerVector isotope_positions(const NumericVector& mz_centroids,
                                const NumericVector& intensity_centroids,
                                const IntegerVector& iso_pattern, double mz,
                                double intensity, double mz_tol,
                                double intensity_tol, int charge, int num_iso) {

    // generate averagine distribution for given mz and charge
    NumericVector intensity_ratios = averagine_ratios(mz, charge, num_iso);

    // initialize vector to collect positions for found isotopes
    IntegerVector positions;

    // iterate over mz_centroids and intensity_centroids to see if they fit
    // expected mz distance and intensity ratios
    int data_size = mz_centroids.size();
    bool keep_looking = true;
    for (int iso = 0; iso < num_iso; ++iso) {
        if (keep_looking) {
            for (int i = 0; i < data_size; ++i) {
                double mz_dist = ppm_diff(mz_centroids[i] - (double)iso / charge, mz);
                double it_dist = frac_diff(intensity_centroids[i]/intensity, intensity_ratios[iso]);
                if (mz_dist < mz_tol && it_dist < intensity_tol && iso_pattern[i] == -1) {
                    positions.push_back(i);
                    break;
                } else if (i == data_size - 1 && mz_dist > mz_tol) {
                    // if no peak for current isotope found then stop looking
                    // for next one since it could be falsly found
                    keep_looking = false;
                }
            }
        }
        else break;
    }
    return positions;
}

//' Convert profile spectra to centroided spectra with assigned charge and
//' position in isotopic envelope.
//'
//' @description \code{charge_singlespec} converts single profile spectrum into
//' centroided spectrum using local maximum across \code{min_width} points.
//' Assigment of charge and isotopic position is based on mz distance between
//' isotopologues and deviation of measured intensity from predicted by
//' averagine model.
//' For high performance, \code{charge_singlespec} is implemented in C++ using
//' \href{http://www.rcpp.org}{Rcpp}.
//'
//' @param spectrum a 2-column numeric matrix with mz [,1] and intenity [,2].
//' @param min_width a size of sliding window for centroiding.
//' @param min_intensity minimal intensity for centroiding.
//' @param mz_tol max mz deviation between isotopologues (ppm).
//' @param intensity_tol max deviation between measured intensity and predicted
//' using averagine model.
//' @param max_charge max charge to assign.
//' @param num_iso max number of isotopologues to search for.
//'
//' @return A \code{data.frame}
//'
//' @export
// [[Rcpp::export]]
DataFrame charge_singlespec(const NumericVector& spectrum, int min_width,
                            double min_intensity, double mz_tol,
                            double intensity_tol, int max_charge, int num_iso) {

  // check arguments
  if ((mz_tol < 0) || (intensity_tol < 0) || (max_charge < 0) || (num_iso < 0)) {
      stop("numerical parameters should be positive");
  }

  // convert profile spectra to centroided
  List centroids = find_centroids(spectrum, min_width, min_intensity);
  NumericVector mz = centroids["mz"];
  NumericVector intensity = centroids["intensity"];

  // initialize vectors for return values
  int data_size = mz.size();
  IntegerVector charge(data_size, 0);
  IntegerVector iso_count(data_size, 1);
  IntegerVector iso_position(data_size, 0);
  IntegerVector iso_pattern(data_size, -1);

  // iterate over centroids and assign charge based on mz distance and intensity ratios
  for (int i = 0; i < data_size; ++i) {
      if (iso_pattern[i] == -1) {
          IntegerVector charge_count(max_charge + 1, 1);
          List charge_positions(max_charge + 1, IntegerVector::create(0));
          for (int charge = max_charge; charge > 0; --charge) {
              IntegerVector current_positions = isotope_positions(mz, intensity, iso_pattern,
                                                                  mz[i], intensity[i], mz_tol, intensity_tol, charge, num_iso);
              charge_count[charge] = current_positions.size();
              charge_positions[charge] = current_positions;
          }
          // find best charge based on number isotopes found
          int best_count = 1;
          int best_charge = 0;
          for (int charge = 0; charge <= max_charge; ++charge) {
              if (charge_count[charge] > best_count) {
                  best_count = charge_count[charge];
                  best_charge = charge;
              }
          }
          // assign patterns to centroids if more that one isotope found
          if (best_count > 1 && best_charge > 0) {
              IntegerVector best_positions = charge_positions[best_charge];
              for (int j = 0; j < best_positions.size(); ++j) {
                  charge[best_positions[j]] = best_charge;
                  iso_count[best_positions[j]] = best_count;
                  iso_position[best_positions[j]] = j;
                  iso_pattern[best_positions[j]] = i;
              }
          }
      }
  }

  return DataFrame::create(_["mz"] = mz,
                           _["intensity"] = intensity,
                           _["charge"] = charge,
                           _["iso_position"] = iso_position,
                           _["iso_pattern"] = iso_pattern,
                           _["iso_count"] = iso_count);
}
