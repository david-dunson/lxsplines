#include <RcppArmadillo.h>
#include <iostream>
using namespace Rcpp;
using namespace std; 

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
NumericVector polyaGamma(arma::mat X, int  truncate){

  int n = X.n_cols;
  arma::mat K(truncate, n, arma::fill::ones);
  arma::mat O = arma::ones(1, K.n_rows);
  double CONST_PI_SQ = 4.0*M_PI*M_PI; 
  double CONST_PI_SQ_HALF = 2.0*M_PI*M_PI;
  for (signed int i = 0; i < truncate;  i++) {
	  K.row(i) += (double) i -0.5; 
  }
  K %= K; 
  
  for (signed int i = 0; i < n; i++){
    K.col(i) +=  X(1,i)*X(1,i)/(CONST_PI_SQ);
  }
  
  K = pow(K, -1.0);

  for (signed int i = 0; i < n; i++) {
	  K.col(i) %= as<arma::vec>(rgamma(truncate, X(0, i), 1));
  }
  
  ////////////////////////////////////////////////////////////////////
  //pgamma = 1.0/pgamma; 
  //pgamma = (1.0/(2.0*M_PI*M_PI))*O*pgamma; 
  //pgamma = pgamma % X.row(0); 
  //arma::mat B = arma::tanh(X.row(1)/2.0)/(2.0*X.row(1));
  //B = B%X.row(0); 
  //B = B/pgamma; 
  /////////////////////////////////////////////////////////////////////

  K = (1.0/(CONST_PI_SQ_HALF))*O*K;
  return wrap(K); 
}
