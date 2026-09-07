#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath> 
#include <limits>
#include <unordered_set>
using namespace Rcpp;

// ====================================================================== //
// ~ Apply two-stream model: used for relating NDVI to reflectance etc. ~ //
// ====================================================================== //
// Function used to solve ndvi for a given pai
double ndvicpp(double pai, double ndviin)
{
  // Red
  double S1 = exp(-0.956782804 * pai);
  double D1 = -0.682265131 / S1 + 0.00796106 * S1;
  double p1 = (0.061666667 / (D1 * S1)) * -0.357497089;
  double p2 = (-0.061666667 * S1 / D1) * 1.556068518;
  double alb_red = p1 + p2;
  // NIR
  double alb_nir = 0.212278228;
  // calculate NDVI
  double ndvi = (alb_nir - alb_red) / (alb_nir + alb_red);
  // Output
  double out = ndvi - ndviin;
  return out;
}
// Root-finding function using the bisection method: ndvi
double solve_ndvi(double ndviin, double tol = 1e-6, int max_iter = 100) {
  double lower = 0.0;
  double upper = 50.0;
  double mid = 0.0;
  for (int iter = 0; iter < max_iter; iter++) {
    mid = (lower + upper) / 2.0;
    double f_lower = ndvicpp(lower, ndviin);
    double f_mid = ndvicpp(mid, ndviin);
    if (std::abs(f_mid) < tol) {
      return mid; // Root found
    }
    if (f_lower * f_mid < 0) {
      upper = mid; // Root lies in the lower half
    }
    else {
      lower = mid; // Root lies in the upper half
    }
  }
  mid = NA_REAL;
  return mid; // Should never reach here
}
// [[Rcpp::export]]
NumericMatrix find_pai(NumericMatrix ndvi)
{
  // Get dimensions
  int rows = ndvi.nrow();
  int cols = ndvi.ncol();
  NumericMatrix pai(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      double val = ndvi(i, j);
      if (!NumericMatrix::is_na(val)) {
        pai(i, j) = solve_ndvi(val);
      }
      else {
        pai(i, j) = NA_REAL;
      }
    }
  }
  return pai;
}
// Function used to calculate leaf or ground reflectance solve
double leafrcpp(double lref, double pai, double gref, double x, double albin)
{
  // Base parameters
  double ltra = 0.33 * lref;
  double om = lref + ltra;
  double a = 1 - om;
  double del = lref - ltra;
  double J = 1.0 / 3.0;
  if (x != 1.0) {
    double mla = 9.65 * pow((3 + x), -1.65);
    if (mla > M_PI / 2) mla = M_PI / 2;
    J = cos(mla) * cos(mla);
  }
  double gma = 0.5 * (om + J * del);
  double h = sqrt(a * a + 2 * a * gma);
  // Calculate base parameters: diffuse
  double S1 = exp(-h * pai);
  double u1 = a + gma * (1 - 1 / gref);
  double D1 = (a + gma + h) * (u1 - h) * 1 / S1 - (a + gma - h)
    * (u1 + h) * S1;
  // Calculate parameters: diffuse
  double p1 = (gma / (D1 * S1)) * (u1 - h);
  double p2 = (-gma * S1 / D1) * (u1 + h);
  double albd = p1 + p2;
  // Output
  double out = albd - albin;
  return out;
}
// Root-finding function using the bisection method: leafr
double solve_lref(double pai, double gref, double x, double albin, double tol = 1e-6, int max_iter = 100) {
  double lower = 0.0001;
  double upper = 0.7499;
  double mid = 0.0;
  for (int iter = 0; iter < max_iter; iter++) {
    mid = (lower + upper) / 2.0;
    double f_lower = leafrcpp(lower, pai, gref, x, albin);
    double f_mid = leafrcpp(mid, pai, gref, x, albin);
    if (std::abs(f_mid) < tol) {
      return mid; // Root found
    }
    if (f_lower * f_mid < 0) {
      upper = mid; // Root lies in the lower half
    }
    else {
      lower = mid; // Root lies in the upper half
    }
  }
  mid = NA_REAL;
  return mid; // Should never reach here
}
// Root-finding function using the bisection method: gref
double solve_gref(double lref, double pai, double x, double albin, double tol = 1e-6, int max_iter = 100) {
  double lower = 0.0001;
  double upper = 0.9999;
  double mid = 0.0;
  for (int iter = 0; iter < max_iter; iter++) {
    mid = (lower + upper) / 2.0;
    double f_lower = leafrcpp(lref, pai, lower, x, albin);
    double f_mid = leafrcpp(lref, pai, mid, x, albin);
    if (std::abs(f_mid) < tol) {
      return mid; // Root found
    }
    if (f_lower * f_mid < 0) {
      upper = mid; // Root lies in the lower half
    }
    else {
      lower = mid; // Root lies in the upper half
    }
  }
  mid = NA_REAL;
  return mid; // Should never reach here
}
// [[Rcpp::export]]
NumericMatrix find_gref(NumericMatrix lref, NumericMatrix pai, 
                        NumericMatrix x, NumericMatrix albin) 
{
  // Get dimensions
  int rows = pai.nrow();
  int cols = pai.ncol();
  NumericMatrix gref(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      double val = pai(i, j);
      if (NumericMatrix::is_na(lref(i, j))) val = NA_REAL;
      if (NumericMatrix::is_na(x(i, j))) val = NA_REAL;
      if (NumericMatrix::is_na(albin(i, j))) val = NA_REAL;
      if (!NumericMatrix::is_na(val)) {
        gref(i, j) = solve_gref(lref(i, j), val, x(i, j), albin(i, j));
      }
      else {
        gref(i, j) = NA_REAL;
      }
    }
  }
  return gref;
}
// [[Rcpp::export]]
NumericMatrix find_lref(NumericMatrix pai, NumericMatrix gref,
                        NumericMatrix x, NumericMatrix albin)
{
  // Get dimensions
  int rows = pai.nrow();
  int cols = pai.ncol();
  NumericMatrix lref(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      double val = pai(i, j);
      if (NumericMatrix::is_na(gref(i, j))) val = NA_REAL;
      if (NumericMatrix::is_na(x(i, j))) val = NA_REAL;
      if (NumericMatrix::is_na(albin(i, j))) val = NA_REAL;
      if (!NumericMatrix::is_na(val)) {
        lref(i, j) = solve_lref(val, gref(i, j), x(i, j),
             albin(i, j));
      }
      else {
        lref(i, j) = NA_REAL;
      }
    }
  }
  return lref;
}
// [[Rcpp::export]]
IntegerMatrix getsoiltypecpp(NumericMatrix bulkden, NumericMatrix clay,
                             NumericMatrix sand, NumericMatrix silt)
{
  double bdensd = 0.08835813;
  double claysd = 0.1378817;
  double sandsd = 0.237984;
  double siltsd = 0.192472;
  NumericVector bulkv = { 1.597779, 1.587082, 1.578984, 1.513506,
                          1.358636, 1.617506, 1.529643, 1.472509, 1.642237, 1.670682,
                          1.604273 };
  NumericVector clayv = { 0, 0.075, 0.15, 0.25, 0.175, 0.275,
                          0.325, 0.275, 0.3, 0.35, 0.5 };
  NumericVector sandv = { 0.5, 0.8, 0.7, 0.5, 0.35, 0.5, 0.75, 0.3, 
                          0.55, 0.1, 0.1 };
  NumericVector siltv = { 0, 0.15, 0.2, 0.3, 0.65, 0.2, 0.3, 0.5, 
                          0.1, 0.1, 0.1 };
  // Loop through and calculate most likely soil type
  // Get dimensions
  int rows = bulkden.nrow();
  int cols = bulkden.ncol();
  IntegerMatrix soiltype(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      double val = bulkden(i, j);
      if (NumericMatrix::is_na(clay(i, j))) val = NA_REAL;
      if (NumericMatrix::is_na(sand(i, j))) val = NA_REAL;
      if (NumericMatrix::is_na(silt(i, j))) val = NA_REAL;
      if (!NumericMatrix::is_na(val)) {
        double difs = 9999.99;
        for (int k = 0; k < 11; ++k) {
          double difn = abs(val - bulkv[k]) / bdensd +
            abs(clay(i, j) - clayv[k]) / claysd +
            abs(sand(i, j) - sandv[k]) / sandsd +
            abs(silt(i, j) - siltv[k]) / sandsd;
          if (difn < difs) soiltype(i, j) = k + 1;
          difs = difn;
          
        } // end k
      } // end NA check
      else {
        soiltype(i, j) = NA_REAL;
      } // end NA check
    } // end j
  } // end i
  return soiltype;
}