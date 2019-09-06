#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <ctime>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

void lmm_2(const arma::vec& y, const arma::mat& W, const arma::mat& X,  const int& maxIter,
              double& sigma2y, double& sigma2beta, arma::vec& beta0, double& loglik_max,
              int& iteration, arma::mat& Sigb, arma::vec& mub){

  int n = y.n_elem, p = X.n_cols;

  if (y.n_elem != X.n_rows || X.n_rows != W.n_rows){
    perror("The dimensions in outcome and covariates (X and W) are not matched");
  }

  if (beta0.n_elem != W.n_cols){
    perror("The dimensions in covariates are not matched in W and beta0");
  }

  if (p != (int)mub.n_elem){
    perror("The dimensions in covariates are not matched in mub");
  }

  if (p != (int)Sigb.n_cols){
    perror("The dimensions in covariates are not matched in Sigb");
  }

  mat XtX = X.t()*X, WtW = W.t()*W, WtX = W.t()*X;
  vec Xty = X.t()*y, Wty = W.t()*y;

  vec SWy;
  mat SWX;

  if(W.n_cols==1){
    SWy = mean(y);
    SWX = mean(X,0);
  } else{
    SWy = solve(WtW, Wty);
    SWX = solve(WtW, WtX);
  }


  double gam, gam2;  // parameter expansion

  vec eVal;
  mat eVec;

  eig_sym(eVal, eVec, XtX);

  // initialize
  sigma2y = var(y);
  sigma2beta = sigma2y/p;
  beta0 = SWy - SWX * mub;
  vec loglik(maxIter);
  loglik(0) = -datum::inf;

  vec D;
  vec Xmu;
  vec y_bar = y - W * beta0;
  double y_Xmu2, E;

  iteration = maxIter-1;
  for (int iter = 1; iter < maxIter; iter ++ ) {
    // E-step
    D = 1 / sigma2beta +  eVal / sigma2y;
    mub = 1/sigma2y * eVec * (eVec.t() * (X.t() * y_bar) / D);
    Xmu = X * mub;
    y_Xmu2 = sum(pow(y_bar-Xmu,2));

    // Evaluate loglik
    E = y_Xmu2/(2*sigma2y) + accu(pow(mub,2))/(2*sigma2beta);
    loglik(iter) = - p*log(sigma2beta)/2 - n*log(sigma2y)/2 - E - sum(log(D))/2 - n/2*log(2*datum::pi);

    if ( loglik(iter) - loglik(iter - 1) < 0 ){
      perror("The likelihood failed to increase!");
    }

    if (abs(loglik(iter) - loglik(iter - 1)) < 1e-10) {
      iteration = iter;
      break;
    }

    // M-step
    gam = sum(y_bar % Xmu) / (accu(pow(Xmu,2)) + sum(eVal/D));
    gam2 = pow(gam , 2);

    beta0 = SWy - (SWX * mub) * gam;
    y_bar = y - W * beta0;;

    sigma2y = sum(pow(y_bar-Xmu*gam,2))/n + gam2 * sum(eVal/D)/n;
    sigma2beta = accu(pow(mub,2))/p + sum(1/D)/p;

    // Reduction step
    sigma2beta = gam2 * sigma2beta;
    // gam = 1;
    // gam2 = pow(gam , 2);
  }

  vec loglik_out;
  loglik_out = loglik.subvec(0, iteration);

  loglik_max = loglik(iteration);
}

// [[Rcpp::export]]
Rcpp::List Estimate_VC_2(const arma::vec y, const arma::mat w, const arma::mat x, const int maxIter){

    double sigma2y = var(y)/2, sigma2beta = var(y)/2, loglik;
    vec beta0 =zeros<vec>(w.n_cols);
    int iter;
    mat Sigb = zeros<mat>(x.n_cols,x.n_cols);
    vec mub  = zeros<vec>(x.n_cols);

    lmm_2(y, w, x, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,Sigb,mub);

    List output = List::create(Rcpp::Named("sigmaerror") = sigma2y,
                               Rcpp::Named("sigmaVC") = sigma2beta,
                               //Rcpp::Named("beta0") = beta0,
                               //Rcpp::Named("loglik_seq") = loglik_out,
                               Rcpp::Named("loglik") = loglik,
                               Rcpp::Named("iteration") = iter,
                               Rcpp::Named("Sigb") = Sigb,
                               Rcpp::Named("mub") = mub);

    return output;

}
