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

void lmm_1(const arma::vec& y, const arma::mat& w, const arma::mat& x, const int& maxIter,
			  double& sigma2y, double& sigma2beta, arma::vec& beta0, double& loglik_max,
			  int& iteration, arma::mat& Sigb, arma::vec& mub){
	if (y.n_elem != x.n_rows || x.n_rows != w.n_rows){
		perror("The dimensions in outcome and covariates (x and w) are not matched");
	}

	sigma2y = var(y), sigma2beta = sigma2y;
    int n = y.n_elem, p = x.n_cols;
	if (beta0.n_elem != w.n_cols){
		perror("The dimensions in covariates are not matched in w and beta0");
	}

	if (p != (int)mub.n_elem){
		perror("The dimensions in covariates are not matched in mub");
	}

	if (p != (int)Sigb.n_cols){
		perror("The dimensions in covariates are not matched in Sigb");
	}

	mat xtx = x.t()*x, wtw = w.t()*w, xtw = x.t()*w;
	vec xty = x.t()*y, wty = w.t()*y;
	double gam, gam2;  // parameter expansion

	vec  dd;
	mat uu;

	eig_sym(dd, uu, xtx);

	mat wtxu = w.t()*x*uu; //q-by-p

	vec dd2;
	mat uu2;

	mat tmp = x*x.t();

	eig_sym(dd2, uu2, tmp);

	vec Cm(p), Cs(p), Cm2(p);
	vec loglik(maxIter);
	vec utxty(p), utxty2(p), utxty2Cm(p), utxty2Cm2(p);
	vec yt(n), yt2(n);

	// evaluate loglik at initial values
	vec uty = uu2.t()*(y - w * beta0);

	vec tmpy = uty / sqrt(dd2 * sigma2beta + sigma2y);
	vec tmpy2 = pow(tmpy,2);
	loglik(0) = -0.5 * sum(log(dd2 * sigma2beta + sigma2y)) - 0.5 * sum(tmpy2);

	//cout << "initial :sigma2beta: " << sigma2beta <<"; sigma2y: " << sigma2y << "; loglik: " << loglik(0) << endl;
	int Iteration = 1;
	for (int iter = 2; iter <= maxIter; iter ++ ) {
		// E-step
		Cm = sigma2y / sigma2beta +  dd;
		Cs = 1 / sigma2beta +  dd / sigma2y;
		Cm2 = pow(Cm , 2);
		// M-step
		utxty = uu.t() * (xty - xtw * beta0); // p-by-1
		utxty2 = pow(utxty, 2);

		utxty2Cm = utxty % utxty / Cm;
		utxty2Cm2 = utxty2Cm/Cm;

		gam = sum(utxty2Cm) / ( sum(dd % utxty2Cm2) + sum(dd / Cs));
		gam2 = pow(gam , 2);
		//cout << " iter: " << iter << "; 1/Cs: " << sum(1 / Cs) << ";utxty2Cm2: " << sum(utxty2Cm2) <<"; gam: " << sum(gam) <<  "; utxty2Cm: " <<	sum(utxty2Cm) <<endl;

		//sigma2beta = ( sum(utxty.t() * diagmat(1 / Cm2) * utxty) + sum(1 / Cs)) / p;
		sigma2beta = ( sum(utxty2Cm2) + sum(1 / Cs)) / p;

		yt = y - w*beta0;
		yt2 = pow(yt , 2);
		sigma2y = (sum(yt2) - 2 * sum(utxty2Cm) * gam + sum(utxty2Cm2 % dd) * gam2 + gam2 * sum(dd / Cs)) / n;

		beta0 = solve(wtw, wty - gam2 * (wtxu * (utxty / Cm)));

		//reduction and reset
		sigma2beta = gam2 * sigma2beta;

		//evaluate loglik
		uty = uu2.t()*(y - w * beta0);
		tmpy = uty / sqrt(dd2 * sigma2beta + sigma2y);
		tmpy2 = pow(tmpy,2);
		loglik(iter - 1) = - 0.5 * sum(log(dd2 * sigma2beta + sigma2y)) - 0.5 * sum(tmpy2);

		//cout << "sigma2beta: " << sigma2beta <<"; sigma2y: " << sigma2y << "; beta0: " << beta0 << "; loglik: " <<	loglik(iter - 1) << "; sum(tmpy2): " << sum(tmpy2) <<endl;
		if ( loglik(iter - 1) - loglik(iter - 2) < 0 ){
			perror("The likelihood failed to increase!");
		}

		Iteration = iter;
		if (abs(loglik(iter - 1) - loglik(iter - 2)) < 1e-10) {

			break;
		}
	}

	Cm = sigma2y / sigma2beta + dd;
	Cs = 1 / sigma2beta + dd / sigma2y;
	Sigb = uu*diagmat(1 / Cs)*uu.t();
	mub = uu * (utxty / Cm);

	vec loglik_out;
	int to = Iteration -1;
	loglik_out = loglik.subvec(0, to);

 	loglik_max = loglik(to);
	iteration = Iteration -1;

}


// [[Rcpp::export]]
Rcpp::List Estimate_VC_1(const arma::vec y, const arma::mat w, const arma::mat x, const int maxIter){

    double sigma2y, sigma2beta, loglik;
    vec beta0 =zeros<vec>(w.n_cols);
    int iter;
    mat Sigb = zeros<mat>(x.n_cols,x.n_cols);
    vec mub  = zeros<vec>(x.n_cols);
    //mat uu;// = zeros<mat>(x.n_cols,x.n_cols);
	//mat D;// = zeros<mat>(x.n_cols,x.n_rows);
    lmm_1(y, w, x, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,Sigb,mub);

    //lmm_pxem_test(y, w, x, maxIter,sigma2y,sigma2beta,beta0,loglik,iter,mub);

    List output = List::create(Rcpp::Named("sigmaerror") = sigma2y,
                               Rcpp::Named("sigmaVC") = sigma2beta,
                               Rcpp::Named("beta0") = beta0,
                               //Rcpp::Named("loglik_seq") = loglik_out,
                               Rcpp::Named("loglik") = loglik,
                               Rcpp::Named("iteration") = iter,
                               Rcpp::Named("Sigb") = Sigb,
                               Rcpp::Named("mub") = mub);
                               //Rcpp::Named("uu") = uu,
							   //Rcpp::Named("D") = D);

    return output;

}
