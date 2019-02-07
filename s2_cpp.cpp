#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double S2_cpp(NumericVector y1, NumericVector y2) {
  
  if (y1.size() != y2.size()) stop("==> y1 and y2 must have the same length.");
  
  long int n = y1.size();
  
  NumericVector diff(pow(n, 2.));
  
  long int k = 0;
  for (register long int i = 0; i < n; ++i){
     for (register long int j = 0; j < n; ++j){
        diff[k] = pow(y1[i] - y2[j], 2.);
        k++;
     }
  } 
  
  NumericVector no_zero = diff[diff != 0];
  
  return Rcpp::median(no_zero);
}
