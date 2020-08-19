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

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


void molmm_pxem_ptr2(const arma::vec& x, const arma::mat& W, const arma::mat& Zx,  const int& maxIter,
                   double& sigma2x, double& sigma2z, arma::vec& mux, double& loglik_max,
                   int& iteration, arma::mat& Sigb, arma::vec& mub){
  
  int n = x.n_elem, p = Zx.n_cols;
  
  if (x.n_elem != Zx.n_rows || Zx.n_rows != W.n_rows){
    perror("The dimensions in outcome and covariates (Zx and W) are not matched");
  }
  
  if (mux.n_elem != W.n_cols){
    perror("The dimensions in covariates are not matched in W and mux");
  }
  
  if (p != (int)mub.n_elem){
    perror("The dimensions in covariates are not matched in mub");
  }
  
  if (p != (int)Sigb.n_cols){
    perror("The dimensions in covariates are not matched in Sigb");
  }
  
  mat ZxtZx = Zx.t()*Zx, WtW = W.t()*W, WtZx = W.t()*Zx;
  vec Zxtx = Zx.t()*x, Wtx = W.t()*x;
  
  vec SWx;
  mat SWZx;
  
  if(W.n_cols==1){
    SWx = mean(x);
    SWZx = mean(Zx,0);
  } else{
    SWx = solve(WtW, Wtx);
    SWZx = solve(WtW, WtZx);
  }
  
  
  double lambda, lambda2;  // parameter expansion
  
  vec eVal;
  mat eVec;
  
  eig_sym(eVal, eVec, ZxtZx);
   
  // initialize
  sigma2x = var(x);
  sigma2z = sigma2x/p;
  mux = SWx - SWZx * mub;
  vec loglik(maxIter);
  loglik(0) = -datum::inf;
  
  vec D;
  vec Zxmu;
  vec x_bar = x - W * mux;
  double x_Zxmu2, E;
  
  iteration = maxIter-1;
  for (int iter = 1; iter < maxIter; iter ++ ) {
    // E-step
    D = 1 / sigma2z +  eVal / sigma2x;
    mub = 1/sigma2x * eVec * (eVec.t() * (Zx.t() * x_bar) / D);
    Zxmu = Zx * mub;
    x_Zxmu2 = sum(pow(x_bar-Zxmu,2));
    
    // Evaluate loglik
    E = x_Zxmu2/(2*sigma2x) + accu(pow(mub,2))/(2*sigma2z);
    loglik(iter) = - p*log(sigma2z)/2 - n*log(sigma2x)/2 - E - sum(log(D))/2 - n/2*log(2*datum::pi);
    
    if ( loglik(iter) - loglik(iter - 1) < 0 ){
      perror("The likelihood failed to increase!");
    }
    
    if (abs(loglik(iter) - loglik(iter - 1)) < 1e-10) {
      iteration = iter;
      break;
    }
    
    // M-step
    lambda = sum(x_bar % Zxmu) / (accu(pow(Zxmu,2)) + sum(eVal/D));
    lambda2 = pow(lambda , 2);
    
    mux = SWx - (SWZx * mub) * lambda;
    x_bar = x - W * mux;
    
    sigma2x = sum(pow(x_bar-Zxmu*lambda,2))/n + lambda2 * sum(eVal/D)/n;
    sigma2z = accu(pow(mub,2))/p + sum(1/D)/p;
    
    // Reduction step
    sigma2z = lambda2 * sigma2z;
    // lambda = 1;
    // lambda2 = pow(lambda , 2);
  }
  
  vec loglik_out;
  loglik_out = loglik.subvec(0, iteration);
  
  loglik_max = loglik(iteration);
}



void mologlike_twas(const arma::mat& R, const arma::vec& res_x, double& res_ysum, const arma::vec& mu, const double& sigma2z, const double& sigma2x,
                  const arma::vec& eVal, const int& n1, const int& n2, const int& p, double& loglik){
  double Eab;//, loglik;
  Eab = sum(res_x % res_x)/2/sigma2x +res_ysum/2 + sum(mu % mu)/2/sigma2z;
  
  loglik = -log(sigma2z)*p*0.5 - log(sigma2x)*n1*0.5 - accu(log(eVal))*n2*0.5 - sum(log(R.diag())) - Eab;
  //cout<<eVal(find(eVal > 0.00001))<<endl;
}



// [[Rcpp::export]]
SEXP moPMR_individual_Egger(SEXP xin, SEXP Yin, SEXP Zxin, SEXP Zyin, SEXP gammain, SEXP alphain, SEXP max_iterin, SEXP epsin){// *
  try{
    const int constrgamma = Rcpp::as<int>(gammain);
    const int constralpha = Rcpp::as<int>(alphain);
    const int maxIter = Rcpp::as<int>(max_iterin);
    const double epsStopLogLik = Rcpp::as<double>(epsin);
    const int pxem_indicator = 1;
    
    const arma::vec x = as<arma::vec>(xin);
    const arma::mat Y = as<arma::mat>(Yin);
    const arma::mat Zx = as<arma::mat>(Zxin);
    const arma::mat Zy = as<arma::mat>(Zyin);
    
    // *convert data format
    int n1 = x.n_elem, n2 = Y.n_cols, p1 = Zx.n_cols, p2 = Zy.n_cols, k = Y.n_rows;
    mat w2;
    mat w1=ones(n1,1);
    if (p1 != p2){
      perror("The dimensions of Zx and Zy are not matched");
    }
    int p = p1;
    if (constrgamma == 1){
      w2 = zeros(1,n2);
    }
    else {
      w2 = sum(Zy.t(),0);
    }

    mat ZxtZx = Zx.t()*Zx, ZytZy = Zy.t()*Zy, w1tw1 = w1.t()*w1, Zxtw1 = Zx.t()*w1;
    vec Zxtx = Zx.t()*x, w1tx = w1.t()*x;
    double w2tw2 = as_scalar(w2*w2.t());
    // initialization using lmm_pxem
    double sigma2x, sigma2z, loglik0;
    vec mux = zeros<vec>(w1.n_cols);
    int iter0;
    mat Sigb = zeros<mat>(Zx.n_cols,Zx.n_cols);
    vec mub  = zeros<vec>(Zx.n_cols);
    
    molmm_pxem_ptr2(x, w1, Zx, maxIter , sigma2x , sigma2z , mux, loglik0 , iter0 , Sigb , mub);
    
    //declaration of variables used within loop
    mat Sigma_inv(p,p), R(p,p), invR(p,p), Sigma(p,p), res_y(k,n2);
    vec mutmp(p), muty = zeros<vec>(p), mu(p), res_x(n1), Zxmu(n1), Zymu(n2),gamma(k);
    mat I = eye<mat>(p,p);

    //initialization of parameters
    if (constrgamma == 1){
      gamma = zeros(k,1);
    }
    else {
      gamma = Y*w2.t()/w2tw2;
    }

    double lambda = 1, lambda2 = lambda*lambda;
    mat alpha = zeros<mat>(k,1), sigmay = cov(Y.t()),S(k,k),sigmay_inv(k,k);

    vec loglik(maxIter), loglik2(maxIter);
    
    // initialization of likelihood
    loglik(0) = NAN;
    int Iteration = 1;
    for (int iter = 2; iter <= maxIter; iter ++ ) {
      // E-step
      vec eVal;
      mat eVec;
      eig_sym(eVal, eVec, sigmay);

      S = diagmat(1.0/eVal);

      sigmay_inv = eVec*S*eVec.t();
      Sigma_inv = ZxtZx/sigma2x + ZytZy*as_scalar(alpha.t()*sigmay_inv*alpha) + I /sigma2z;

      R = chol(Sigma_inv);
      invR = inv(R);
      Sigma = invR*invR.t();

      mutmp = (Zxtx - Zxtw1*mux)/sigma2x + Zy.t()*(Y -gamma*w2).t()*sigmay_inv*alpha;
      mu = Sigma*mutmp;

      //evaluate incomplete log-likelihood
      Zxmu = Zx*mu;
      Zymu = Zy*mu;
      res_x = x - Zxmu*lambda - w1*mux;
      res_y = Y - alpha*Zymu.t() -gamma*w2;
      
      mat s= sigmay_inv*res_y*res_y.t();
      double res_ysum = sum(s.diag());
      //loglik(iter - 1) = loglike_twas(R, res_x, res_y, mu, sigma2z, sigma2x, sigma2y,n1,n2,p);
      mologlike_twas(R, res_x, res_ysum, mu, sigma2z, sigma2x, eVal, n1 , n2, p, loglik(iter - 1));
      
      if ( loglik(iter - 1) - loglik(iter - 2) < -1e-7){
        perror("The likelihood failed to increase!");
      }
      // M-step
      mat trtmp1 = ZxtZx * Sigma;
      double tr1 = trace(trtmp1);
      
      if (pxem_indicator == 1){
        lambda = sum( (Zxtx - Zxtw1*mux) % mu )/(sum(Zxmu % Zxmu) + tr1);
        lambda2 = lambda*lambda;
      }
      else {
        lambda = 1; lambda2 = 1;
      }
      
      if (constrgamma == 1){
         gamma = zeros(k,1);
      }
      else {
        gamma = (Y-alpha*(Zymu.t()))*w2.t()/w2tw2;
      }
      mat trtmp2 = ZytZy * Sigma;
      double tr2 = trace(trtmp2);
      if (constralpha == 1){
        mat alpha = zeros<mat>(k,1);
      }
      else {
        alpha = (Y -gamma*w2)* Zymu/(sum(Zymu % Zymu) + tr2);
      }
      mat alpha2(k,k);
      alpha2 = alpha*(alpha.t());
      
      res_x = x - Zxmu*lambda - w1*mux;
      res_y = Y - alpha*(Zymu.t()) -gamma*w2;
      
      sigma2x = (sum(res_x % res_x) + lambda2 * tr1)/n1;
      sigmay = ((res_y * res_y.t()) + alpha2 * tr2)/n2;
      sigma2z = (sum(mu % mu) + trace(Sigma))/p;

      // Reduction-step
      sigma2z = lambda2 * sigma2z;
      alpha = alpha / lambda;
      alpha2 =  alpha*alpha.t();
      lambda = 1;
      lambda2 = 1;
      
      Iteration = iter;
      if (iter > 2){
        if (abs(loglik(iter - 1) - loglik(iter - 2)) < epsStopLogLik) {
          
          break;
        }
      }
    }
    
    
    vec loglik_out;
    int to = Iteration -1;
    loglik_out = loglik.subvec(0, to);
    
    double loglik_max = loglik(to);
    
    
    if (p * sigma2z < 1e-04){
      perror("The estimate of gene expression heritability explained by cis-SNPs is smaller than 0.01%");
    }
    
    
    
    return List::create(Rcpp::Named("alpha") = alpha,
                        Rcpp::Named("gamma") = gamma,
                        Rcpp::Named("sigmaX") = sigma2x,
                        Rcpp::Named("sigmaY") = sigmay,
                        Rcpp::Named("sigmaZ") = sigma2z,
                        Rcpp::Named("loglik_seq") = loglik_out,
                        Rcpp::Named("loglik") = loglik_max,
                        Rcpp::Named("iteration") = Iteration-1);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
  }
  return R_NilValue;
}// end func
