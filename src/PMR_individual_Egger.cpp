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


void lmm_pxem_ptr2(const arma::vec& y, const arma::mat& W, const arma::mat& X,  const int& maxIter,
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



void loglike_twas(const arma::mat& R, const arma::vec& res_y, const arma::vec& res_z, const arma::vec& mu, const double& sigma2beta, const double& sigma2y,
					const double& sigma2z, const int& n1, const int& n2, const int& p, double& loglik){
	double Eab;//, loglik;
	Eab = sum(res_y % res_y)/2/sigma2y + sum(res_z % res_z)/2/sigma2z + sum(mu % mu)/2/sigma2beta;

	loglik = -log(sigma2beta)*p*0.5 - log(sigma2y)*n1*0.5 - log(sigma2z)*n2 *0.5 - sum(log(R.diag())) - Eab;

}



// [[Rcpp::export]]
SEXP PMR_individual_Egger(SEXP yin, SEXP zin, SEXP x1in, SEXP x2in, SEXP gammain, SEXP alphain, SEXP max_iterin, SEXP epsin){// *
try{
	const int constrgamma = Rcpp::as<int>(gammain);
	const int constralpha = Rcpp::as<int>(alphain);
	const int maxIter = Rcpp::as<int>(max_iterin);
	const double epsStopLogLik = Rcpp::as<double>(epsin);
	const int pxem_indicator = 1;
	
	const arma::vec y = as<arma::vec>(yin);
	const arma::vec z = as<arma::vec>(zin);
	const arma::mat x1 = as<arma::mat>(x1in);
	const arma::mat x2 = as<arma::mat>(x2in);
	
	// *convert data format
	int n1 = y.n_elem, n2 = z.n_elem, p1 = x1.n_cols, p2 = x2.n_cols;
    mat w2;
    mat w1=ones(n1,1);
    if (p1 != p2){
        perror("The dimensions of x1 and x2 are not matched");
    }
    int p = p1;
       if (constrgamma == 1){
       // w2= ones(n2,1);
	   w2= zeros<mat>(n2,1);
   }
  else {
            w2 = sum(x2,1);
             }
    mat x1tx1 = x1.t()*x1, x2tx2 = x2.t()*x2, w1tw1 = w1.t()*w1, w2tw2 = w2.t()*w2, x1tw1 = x1.t()*w1, x2tw2 = x2.t()*w2;
    vec x1ty = x1.t()*y, x2tz = x2.t()*z, w2tz = w2.t()*z, w1ty = w1.t()*y;

    // initialization using lmm_pxem
    double sigma2y, sigma2beta, loglik0;
    vec beta0 =zeros<vec>(w1.n_cols);
    int iter0;
    mat Sigb = zeros<mat>(x1.n_cols,x1.n_cols);
    vec mub  = zeros<vec>(x1.n_cols);

    lmm_pxem_ptr2(y, w1, x1, maxIter,sigma2y,sigma2beta,beta0,loglik0,iter0,Sigb,mub);

    //declaration of variables used within loop
    mat Sigma_inv(p,p), R(p,p), invR(p,p), Sigma(p,p);
    vec mutmp(p), mu(p), res_y(n1), res_z(n2), x1mu(n1), x2mu(n2);
    double ztmp;
    mat I = eye<mat>(p,p);
    double alpha0;
    //initialization of parameters
	if (constrgamma == 1){
       // w2= ones(n2,1);
	    alpha0 =0;
   }
   else {
	    alpha0 = as_scalar(solve(w2tw2, w2tz));           
             }
   // vec alpha0 = solve(w2tw2, w2tz);
    double alpha = 0, sigma2z = var(z), gam = 1, gam2 = gam*gam, alpha2 = alpha*alpha;

    vec loglik(maxIter), loglik2(maxIter);

    // initialization of likelihood
    loglik(0) = NAN;
    int Iteration = 1;
    for (int iter = 2; iter <= maxIter; iter ++ ) {
        // E-step
        Sigma_inv = x1tx1/sigma2y + x2tx2*(alpha2/sigma2z) + I /sigma2beta;
        //cout << "break : " << Sigma_inv << endl;
        R = chol(Sigma_inv);
        invR = inv(R);
        Sigma = invR*invR.t();
        mutmp = (x1ty - x1tw1*beta0)/sigma2y + (x2tz - x2tw2*alpha0)*(alpha/sigma2z);
        mu = Sigma*mutmp;

        //evaluate incomplete log-likelihood
        x1mu = x1*mu;
        x2mu = x2*mu;
        res_y = y - x1mu*gam - w1*beta0;
        res_z = z - alpha*x2mu - w2*alpha0;

        //loglik(iter - 1) = loglike_twas(R, res_y, res_z, mu, sigma2beta, sigma2y, sigma2z,n1,n2,p);
        loglike_twas(R, res_y, res_z, mu, sigma2beta, sigma2y, sigma2z,n1,n2,p, loglik(iter - 1));

        if ( loglik(iter - 1) - loglik(iter - 2) < -1e-7){
            perror("The likelihood failed to increase!");
        }
        // M-step
        mat trtmp1 = x1tx1 * Sigma;
        double tr1 = trace(trtmp1);

        if (pxem_indicator == 1){
            gam = sum( (x1ty - x1tw1*beta0) % mu )/(sum(x1mu % x1mu) + tr1);
            gam2 = gam*gam;
        }
        else {
            gam = 1; gam2 = 1;
        }

        beta0 = solve(w1tw1, w1ty - gam*x1tw1.t()*mu);
		if (constrgamma == 1){
       // w2= ones(n2,1);
	    alpha0 =0;
   }
   else {
	   alpha0= as_scalar(solve(w2tw2, w2tz - alpha*(x2tw2.t()*mu)));           
             }
		
        //alpha0= solve(w2tw2, w2tz - alpha*(x2tw2.t()*mu));

        ztmp = sum( (z - w2*alpha0)% x2mu);
        mat trtmp2 = x2tx2 * Sigma;
        double tr2 = trace(trtmp2);
        if (constralpha == 1){
            alpha = 0;
        }
        else {
            alpha = ztmp/(sum(x2mu % x2mu) + tr2);
        }
        alpha2 = alpha*alpha;

        res_y = y - x1mu*gam - w1*beta0;
        res_z = z - alpha*x2mu - w2*alpha0;

        sigma2y = (sum(res_y % res_y) + gam2 * tr1)/n1;
        sigma2z = (sum(res_z % res_z) + alpha2 * tr2)/n2;
        sigma2beta = (sum(mu % mu) + trace(Sigma))/p;

        // Reduction-step
        sigma2beta = gam2 * sigma2beta;
        alpha = alpha / gam;
        alpha2 = alpha*alpha;
        gam = 1;
        gam2 = 1;

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
    
   
	  if (p * sigma2beta < 1e-04){
    perror("The estimate of gene expression heritability explained by cis-SNPs is smaller than 0.01%");
  }
  
 
	
    return List::create(Rcpp::Named("alpha") = alpha,
                               Rcpp::Named("gamma") = alpha0,
                               Rcpp::Named("sigmaX") = sigma2y,
                                Rcpp::Named("sigmaY") = sigma2z,
                               Rcpp::Named("sigmabeta") = sigma2beta,
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
