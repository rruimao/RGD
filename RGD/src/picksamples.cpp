//#include<Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;


// [[Rcpp::export]]
Rcpp::List picksamples(const arma::mat X,const arma::mat Y,const double n_ratio){
  int n=X.n_rows;
  int p=X.n_cols-1;
  int n_sp=round(n*n_ratio);
  arma::mat X_sp(n_sp, p+1, fill::zeros);
  arma::mat Y_sp(n_sp,1,fill::zeros);
  uvec index=randperm(n,n_sp);
  for (int i=0;i<n_sp;i++){
    int j=index[i];
    Y_sp[i]=Y[j];
    X_sp.row(i)=X.row(j);
  }
return Rcpp::List::create(Rcpp::Named("X_sp")=X_sp,Rcpp::Named("Y_sp")=Y_sp);
  }

