//#include<Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;


// [[Rcpp::export]]
arma::vec GD(arma::vec beta_0,const double tau,const arma::mat X,const arma::mat Y,double eta_0=1e-3,const double alpha=2){
  int n=X.n_rows;
  arma::vec e_0=Y-X*beta_0;
  arma::vec a=arma::abs(e_0);
  uvec b=arma::find(a>tau);
  double L_0=0.5*arma::sum(a.elem(arma::find(a<tau))%a.elem(arma::find(a<tau)))+0.5*tau*tau*b.n_rows;
  arma::mat w(n, 1, fill::ones);
  w.elem(arma::find(e_0<-tau)).zeros();
  w.elem(arma::find(e_0>tau)).zeros();
  arma::vec ew=e_0%w;
  arma::vec diff=-X.t()*ew/n;
  arma::vec beta=beta_0-eta_0*diff;
  arma::vec e_1=Y-X*beta;
  a=arma::abs(e_1);
  b=arma::find(a>tau);
  double L_1=0.5*arma::sum(a.elem(arma::find(a<tau))%a.elem(arma::find(a<tau)))+0.5*tau*tau*b.n_rows;
  double eta=alpha*eta_0;
  beta=beta_0-eta*diff;
  arma::vec e_2=Y-X*beta;
  a=arma::abs(e_2);
  b=arma::find(a>tau);
  double L_2=0.5*arma::sum(a.elem(arma::find(a<tau))%a.elem(arma::find(a<tau)))+0.5*tau*tau*b.n_rows;
  
  while (L_2<L_1){
    L_1=L_2;
    eta=alpha*eta;
    beta=beta_0-eta*diff;
    e_2=Y-X*beta;
    a=arma::abs(e_2);
    b=arma::find(a>tau);
    L_2=0.5*arma::sum(a.elem(arma::find(a<tau))%a.elem(arma::find(a<tau)))+0.5*tau*tau*b.n_rows;
  }
  eta=eta/alpha;
  beta=beta_0-eta*diff;
  if (L_0<L_1){
    beta=beta_0;
  }

  int j=0;
  double d= sqrt(arma::sum((beta_0-beta)%(beta_0-beta))/n);
  while (d > 1e-8  &&  j<10000000){
    beta_0=beta;
    e_0=Y-X*beta_0;
    a=arma::abs(e_0);
    b=arma::find(a>tau);
    L_0=0.5*arma::sum(a.elem(arma::find(a<tau))%a.elem(arma::find(a<tau)))+0.5*tau*tau*b.n_rows;
    arma::mat w(n, 1, fill::ones);
    w.elem(arma::find(e_0<-tau)).zeros();
    w.elem(arma::find(e_0>tau)).zeros();
    ew=e_0%w;
    diff=-X.t()*ew/n;
    beta=beta_0-eta_0*diff; 
    e_1=Y-X*beta;
    a=arma::abs(e_1);
    b=arma::find(a>tau);
    L_1=0.5*arma::sum(a.elem(arma::find(a<tau))%a.elem(arma::find(a<tau)))+0.5*tau*tau*b.n_rows;
    eta=alpha*eta_0;
    beta=beta_0-eta*diff;
    e_2=Y-X*beta;
    a=arma::abs(e_2);
    b=arma::find(a>tau);
    L_2=0.5*arma::sum(a.elem(arma::find(a<tau))%a.elem(arma::find(a<tau)))+0.5*tau*tau*b.n_rows;
    while (L_2<L_1){
      L_1=L_2;
      eta=alpha*eta;
      beta=beta_0-eta*diff;
      e_2=Y-X*beta;
      a=arma::abs(e_2);
      b=arma::find(a>tau);
      L_2=0.5*arma::sum(a.elem(arma::find(a<tau))%a.elem(arma::find(a<tau)))+0.5*tau*tau*b.n_rows;
    }
    eta=eta/alpha;
    beta=beta_0-eta*diff;
    if (L_0<L_1){
      beta=beta_0;
    }
    d= sqrt(arma::sum((beta_0-beta)%(beta_0-beta))/n);
    j++;
  }
  return beta;
}


// [[Rcpp::export]]
arma::mat runif_in_pball(const int n,const int d, const int p, const double r=1){
  arma:: mat out(n, d, fill::ones);
  const int m=n*d;
  double l=floor(1.0/p);
  double p1=1.0/p-l;
  double p2=1-p1;
  int i=0;
  arma::vec R(m,fill::zeros);
  for (;i<m;++i){
    arma::vec ss(l,fill::randu);  
    double s=arma::sum(log(ss));
    if (l!=1.0/p){
      arma::vec U(2, fill::randu);
      double a=pow(U[1],(1.0/p1))+pow(U[2],(1.0/p2));
      while (a>1){
        arma::vec U(2, fill::randu);
      }
      arma::vec A(1,fill::randu);
      s=s+sum(log(A))*pow(U[0],1.0/p1)/(pow(U[0],(1.0/p1))+pow(U[1],(1.0/p2)));
    }
    R[i]=pow((-p*s),1.0/p);
  }
  arma::vec A = randi<arma::vec>(m, distr_param(0, 1));
  arma::vec B=A-1;
  arma::vec C=A+B;
  arma::mat epsilon=R%C;
  epsilon.set_size(n,d);
  arma::mat signs=2*randi<arma::vec>(m, distr_param(1, 2))-1;
  signs.set_size(n, d);
  arma::mat x=signs%epsilon;
  arma::vec z=pow(randu(n),1.0/d);
  for (int j=0;j<n;++j){
    out.row(j)=r*z[j]*x.row(j)/pow(accu(abs(pow(x.row(j),p))),(1.0/p));
  }
  return out;
}



// [[Rcpp::export]]
arma::vec RGD(const arma::mat X,const arma::mat Y,const double tau,const double iter,double eta_0=1e-3,const double alpha=2){
  arma::vec beta_tmp=(X.t()*X).i()*X.t()*Y;
  arma::vec beta_0=(X.t()*X).i()*X.t()*Y;
  arma::vec beta=(X.t()*X).i()*X.t()*Y;
  int p=X.n_cols-1;
  double L=1000000000;
  for (int s=0; s<iter;s++){
    beta_0=runif_in_pball(1,p+1,2,tau).t();
    arma::vec beta_tmp=GD(beta_0,tau,X,Y,eta_0,alpha);
    arma::vec r=Y-X*beta_tmp;
    arma::vec a=arma::abs(r);
    uvec b=arma::find(a>tau);
    double L_tmp=0.5*arma::sum(a.elem(arma::find(a<tau))%a.elem(arma::find(a<tau)))+0.5*tau*tau*b.n_rows;
    if (L_tmp<L){
      L=L_tmp;
      beta=beta_tmp;
    }
  }
  return beta;
}


