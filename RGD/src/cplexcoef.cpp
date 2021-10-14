//#include<Rcpp.h>
#include <RcppArmadillo.h>
#include <utility>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace std;


// [[Rcpp::export]]
Rcpp::List cplexcoef(const arma::mat X,const arma::mat Y,const double tau){
  int n=X.n_rows;
  int p=X.n_cols-1;
  arma::mat Q(3*n+p+1, 3*n+p+1, fill::zeros);
  for (int i=0;i<n;i++){
    Q(i,i)=2;
  }
  arma::mat f(1,3*n+p+1,fill::zeros);
  for (int j=2*n;j<3*n;j++){
    f.col(j)=pow(tau,2);
  }
 int  M=10000;
 arma::mat I(n,n); 
 I.eye();
  arma::mat Z(n,n,fill::zeros);
  arma::mat Z1(n,p+1,fill::zeros);
  arma::mat A1=join_rows(-I,-I,Z,-X);
  arma::mat A2=join_rows(-I,-I,Z,X);
  arma::mat A3=join_rows(Z,I,-M*I,Z1);
  arma::mat A=join_cols(A1,A2,A3);
  arma::mat Z2(n,1,fill::zeros);
  arma::mat b=join_cols(-Y,Y,Z2);
  arma::mat O(n+p+1,1,fill::ones);
  arma::mat O1(n,1,fill::ones);
  arma::mat O2(2*n+p+1,1,fill::ones);
  arma::mat lb=join_cols(Z2,Z2,-datum::inf*O);
  arma::mat ub=join_cols(tau*O1,datum::inf*O2);
// 1 for continuous and 2 for binary
  arma::mat vtype(3*n+p+1,1,fill::zeros);
  for (int i=0;i<2*n;++i){
    vtype[i]=1;
  }
  for (int i=2*n;i<3*n;++i){
    vtype[i]=2;
  }
  for (int i=3*n;i<3*n+p+1;++i){
    vtype[i]=1;
  }
return Rcpp::List::create(Rcpp::List::create(Rcpp::Named("Q")=Q,Rcpp::Named("f")=f,
                            Rcpp::Named("A")=A,Rcpp::Named("b")=b,
                            Rcpp::Named("lb")=lb,Rcpp::Named("ub")=ub,
                            Rcpp::Named("vtype")=vtype));
  
}


