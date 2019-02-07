#include <Rcpp.h>
#include <math.h>
#include <omp.h>

// [[Rcpp::plugins(openmp)]]

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
double S2(NumericVector y1, NumericVector y2) {
  
  if (y1.size() != y2.size()) stop("==> y1 and y2 must have the same length.");
  
  int n = y1.size();
  
  NumericVector diff(pow(n, 2.));
  
  int k = 0;
  
  for (int i = 0; i < n; ++i){
     for (int j = 0; j < n; ++j){
        diff[k] = pow(y1[i] - y2[j], 2.);
        k++;
     }
  } 
  
  NumericVector no_zero = diff[diff != 0];
  
  return Rcpp::median(no_zero);
  
  //return median(y1);
}

/*** R
S2R = function(y1, y2){
   D = outer(y1, y2, '-')
   D = D^2
   D_no_zero = D[which(!D == 0)]
   median(D_no_zero)
}

# y1 <- y2 <- 1:3000
# 
#microbenchmark::microbenchmark(S2(y1,y2))
#microbenchmark::microbenchmark(S2R(y1,y2))

#S2(y1,y2)
#S2R(y1,y2)

*/

