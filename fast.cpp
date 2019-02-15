# Author: Pedro Rafael D. Marinho
# Institution: Federal University of Paraíba, Brazil
# Department of Statistics
# Email: pedro@de.ufpb.br


//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/


#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]
double S2_cpp(NumericVector y1, NumericVector y2) {
  
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

// [[Rcpp::export]]
double det_cpp(arma::mat& x){
   double result = arma::det(x);
   return result;
}

// [[Rcpp::export]]
NumericVector diag_cpp(arma::mat& x){
   arma::vec result = arma::diagvec(x);
   return NumericVector(result.begin(),result.end());
}

// [[Rcpp::export]]
arma::mat solve_cpp(arma::mat& x){
   double determinant = arma::det(x);
   
   arma::mat result;
   
   if (determinant != 0){
      result = arma::inv(x);
   }else{
      /* Moore-Penrose pseudo-inverse of matrix x */
      result = arma::pinv(x);
   }
   
   return result;
}

// /*** R
// S2R = function(y1, y2){
//    D = outer(y1, y2, '-')
//    D = D^2
//    D_no_zero = D[which(!D == 0)]
//    median(D_no_zero)
// }
// 
//  y1 <- y2 <- 1:1e4
// 
// #microbenchmark::microbenchmark(S2_cpp(y1,y2), times = 5L)
// #microbenchmark::microbenchmark(S2R(y1,y2), times = 5L)
// 
// S2_cpp(y1,y2) # Aqui é C++.
// S2R(y1,y2) # Aqui é R.

