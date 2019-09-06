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

void lmm_pxem_ptr2_summary(const arma::vec& betaxh, const arma::mat& Sigma1,  const int& n1, const int& maxIter,
                   double& sigma2x, double& sigma2beta, double& loglik_max,
                   int& iteration, arma::mat& Sigb, arma::vec& mub){
  
  int n = betaxh.n_elem, p = Sigma1.n_cols;
  //total number of elements; present in Mat, Col, Row, Cube, field and SpMat
  //number of rows; present in Mat, Col, Row, Cube, field and SpMat
  if (n != p){
    perror("The dimensions in betaxh and Sigma1 are not matched");//Print error message
  }

   if (p != (int)mub.n_elem){
    perror("The dimensions in covariates are not matched in mub");
  }
  
  if (p != (int)Sigb.n_cols){
    perror("The dimensions in covariates are not matched in Sigb");
  }
  

  double lambda, lambda2;  // parameter expansion

 
  
  vec eVal;
  mat eVec;
  
  eig_sym(eVal, eVec, Sigma1);
  //sigma2x = var(betaxh)/ sum(1/eVal);
  sigma2x = 1;
  sigma2beta = sigma2x/p;

  vec loglik(maxIter);
  loglik(0) = -datum::inf;
  //loglik:最大似然  datum::inf:∞, infinity
  vec D;
  double y_Xmu2,E;
  
  iteration = maxIter-1;
  for (int iter = 1; iter < maxIter; iter ++ ) {
    // E-step
    D = 1 / sigma2beta +  eVal / sigma2x;
	
    // mub = 1/sigma2x * eVec * (eVec.t()* Sigma1*(betaxh-mub) / D);
    
	mub = 1/sigma2x * eVec * (eVec.t()* (n1-1)* betaxh / D);
	
	y_Xmu2=as_scalar((n1-1)-2*(n1-1)*mub.t()*betaxh + mub.t()*Sigma1*mub);
	
    // Evaluate loglik
	
    E = as_scalar(y_Xmu2/(2*sigma2x) + accu(pow(mub,2))/(2*sigma2beta));
	
    loglik(iter) = - p*log(sigma2beta)/2 - n1*log(sigma2x)/2 - E - sum(log(D))/2 - n1/2*log(2*datum::pi);
	
	//loglik(iter) = - p*log(sigma2beta)/2 - n1*log(sigma2x)/2 - E - sum(log(D))/2;
    
   //accu()Accumulate (sum) all elements of a vector, matrix or cube
    //datum::pi	π, the ratio of any circle's circumference to its diameter
    if ( loglik(iter) - loglik(iter - 1) < 0 ){
      perror("The likelihood failed to increase!");
	  //cout<<"here2"<<endl;
    }
    
    if (abs(loglik(iter) - loglik(iter - 1)) < 1e-10) {
      iteration = iter;
      break;
    }
    
    // M-step
    lambda = as_scalar( (mub.t()*(n1-1)*betaxh)/ (mub.t()*Sigma1*mub + sum(eVal/D)));
	
    lambda2 = pow(lambda , 2);
    //%	Schur product: element-wise multiplication of two objects
	//vec res_betaxh;
    //res_betaxh = betaxh-lambda*mub;
	double res_betaxh;
	
	res_betaxh=as_scalar((n1-1)-2*(n1-1)*lambda*mub.t()*betaxh+lambda2*mub.t()*Sigma1*mub);
	
	sigma2x = as_scalar((res_betaxh + lambda2 * sum(eVal/D))/n1);
	
    //sigma2x = as_scalar(((res_betaxh.t()*Sigma1*res_betaxh) + lambda2 * sum(eVal/D))/p);
	
	
    //sigma2beta = accu(pow(mub,2))/p + sum(Sigma)/p;
    sigma2beta = accu(pow(mub,2))/p + sum(1/D)/p;
    // Reduction step
    sigma2beta = lambda2 * sigma2beta;
    // gam = 1;
    // gam2 = pow(gam , 2);
  }
  
  vec loglik_out;
  loglik_out = loglik.subvec(0, iteration);
  
  loglik_max = loglik(iteration);
}

void loglike_twas_summary(const arma::mat& R, const arma::mat& Sigma1, const arma::mat& Sigma2, const arma::vec& res_betaxh, const arma::vec& res_betayh,const arma::vec& mu,
                  const double& sigma2beta, const double& sigma2x, const double& sigma2y, const int& n1, const int& n2, const int& p, double& loglik){
  double Eab;//, loglik;
  // Eab = sum(res_betayh.t()* Sigma2 * res_betayh)/2/sigma2y + sum(res_betaxh.t() * Sigma1 * res_betaxh)/2/sigma2x + sum(mu % mu)/2/sigma2beta;
  
  Eab = as_scalar(res_betayh/2/sigma2y) + as_scalar(res_betaxh)/2/sigma2x + sum(mu % mu)/2/sigma2beta;
  
  loglik = -log(sigma2beta)*p*0.5 - log(sigma2x)*n1*0.5 - log(sigma2y)*n2 *0.5 - sum(log(R.diag())) - Eab;
  
}




// [[Rcpp::export]]

SEXP PMR_summary_Egger_CPP(SEXP betaxin, SEXP betayin, SEXP Sigma1sin, SEXP Sigma2sin, SEXP samplen1, SEXP samplen2, SEXP gammain, SEXP alphain, SEXP max_iterin, SEXP epsin){// *
  try{ 
    const int constrgamma = Rcpp::as<int>(gammain);
    const int constralpha = Rcpp::as<int>(alphain);
    const int n1 = Rcpp::as<int>(samplen1);  
    const int n2 = Rcpp::as<int>(samplen2); 
    const int maxIter = Rcpp::as<int>(max_iterin);
    const double epsStopLogLik = Rcpp::as<double>(epsin);
    const int pxem_indicator = 1;
    		
    const arma::vec betax = as<arma::vec>(betaxin);
    const arma::vec betay = as<arma::vec>(betayin);
    const arma::mat Sigma1s = as<arma::mat>(Sigma1sin);
    const arma::mat Sigma2s = as<arma::mat>(Sigma2sin);
    
    // *convert data format
    int p1 = Sigma1s.n_cols, p2 = Sigma2s.n_cols;
   mat w2;
   //vec w2;
    if (p1 != p2){
      perror("The dimensions of Sigma1 and Sigma2 are not matched");
    }
    int p = p1;
    if (constrgamma == 1){
       w2= zeros<mat>(p2,1);
    }
    else {
      w2 = ones(p2,1);
    }
	
	vec betaxh=betax, betayh=betay;
	  
  mat Sigma1 = Sigma1s* (n1-1);
  mat Sigma2 = Sigma2s* (n2-1);
  //mat Sigma1 = Sigma1s;
  //mat Sigma2 = Sigma2s;
  // mat Sigma1 = Sigma1s;
  // mat Sigma2 = Sigma2s;
   
  // vec Sigma2tbetayh=Sigma2.t()*betayh, Sigma2tw2=Sigma2.t()*w2,Sigma1tbetaxh=Sigma1.t()*betaxh,w2tbetayh=w2.t()*betayh;
  
    // initialization using lmm_pxem
    double sigma2x, sigma2beta, loglik0;
    int iter0;
    mat Sigb = zeros<mat>(p,p);
    vec mub  = zeros<vec>(p);
    //sigma2x = 0.9;
	//sigma2beta = 0.1/p;
     lmm_pxem_ptr2_summary(betaxh, Sigma1, n1, maxIter,sigma2x,sigma2beta,loglik0,iter0,Sigb,mub);
    
    
    mat Sigma_inv(p,p), R(p,p), invR(p,p), Sigma(p,p);
    vec mutmp(p), mu(p), res_betaxh(p), res_betayh(p), Sigma1mu(p), Sigma2mu(p);
    double ztmp;
    mat I = eye<mat>(p,p);
    // eye<mat>(p,p)
    //initialization of parameters
    // vec gamma = mean(betayh);
	double gamma;
	
	if (constrgamma == 1){
       // w2= ones(n2,1);
	    gamma =0;
   }
   else {
	    //gamma = as_scalar(mean(betayh));  
         mat ww;
		 ww=ones<mat>(1,p2);
		 gamma=as_scalar((ww * betayh * (n2-1))/accu(Sigma2));		
             }
	
//	double gamma = mean(betayh);
	
	//cout<< gamma << endl;
	 
	vec eVal;
    mat eVec;
    eig_sym(eVal, eVec, Sigma2);
    //double alpha = 0, sigma2y =  var(betayh)/ sum(1/eVal), lambda = 1, lambda2 = lambda*lambda, alpha2 = alpha*alpha;
   double alpha = 0, sigma2y =  1, lambda = 1, lambda2 = lambda*lambda, alpha2 = alpha*alpha;
   
    vec loglik(maxIter), loglik2(maxIter);
    
    // initialization of likelihood
    loglik(0) = NAN;
    int Iteration = 1;
    for (int iter = 2; iter <= maxIter; iter ++ ) {
      // E-step
      Sigma_inv = Sigma1/sigma2x + Sigma2*(alpha2/sigma2y) + I /sigma2beta;
      //cout << "break : " << Sigma_inv << endl;
      R = chol(Sigma_inv);
      invR = inv(R);
      Sigma = invR*invR.t();
	 // cout<< Sigma <<endl;
      mutmp = (n1-1)*betaxh/sigma2x + ((n2-1)*betayh - sum(Sigma2,1)*gamma)*(alpha/sigma2y);
	  // cout<<mutmp<<endl;
      mu = Sigma*mutmp;
      
      //evaluate incomplete log-likelihood
      //res_betaxh = betaxh - mu*lambda;
	  
	  res_betaxh=as_scalar((n1-1)-2*(n1-1)*lambda*mu.t()*betaxh+lambda2*mu.t()*Sigma1*mu);
	    
	 //cout<<res_betaxh<<endl;	
		
      res_betayh = as_scalar((n2-1)-2*alpha*(n2-1)*mu.t()*betayh-2*gamma*(n2-1)*sum(betayh)+ alpha2*mu.t()*Sigma2*mu + gamma*gamma*accu(Sigma2)+2*alpha*gamma*mu.t()*sum(Sigma2,1));
	  
     
     loglike_twas_summary(R,Sigma1,Sigma2,res_betaxh, res_betayh,mu,sigma2beta, sigma2x, sigma2y,n1,n2,p, loglik(iter - 1));
      //lambda, alpha,mu,Sigma1,Sigma2, Sigma,sigma2beta, sigma2x, sigma2y,p,
	 // cout<<"after loglik i-1 = " << loglik(iter - 1) << endl;
      if ( loglik(iter - 1) - loglik(iter - 2) < -1e-7){
        perror("The likelihood failed to increase!");
		 //cout<<"here"<<endl;
      }
      // M-step
      mat trtmp1 = Sigma1 * Sigma;
      double tr1 = trace(trtmp1);
      
      if (pxem_indicator == 1){
        lambda = as_scalar(( betaxh.t() *(n1-1)* mu )/((mu.t() * Sigma1 * mu) + tr1));
        lambda2 = lambda*lambda;
      }
      else {
        lambda = 1; lambda2 = 1;
      }
//      gamma= solve(w2tw2, w2tbetayh - alpha*mu);
      // if (w2.t()*w2 == 0){
	if (constrgamma == 1){
       gamma= 0;
      }
       else {
         // gamma= sum(sum(Sigma2,0)%(betayh - alpha*mu))/accu(Sigma2);
		 // gamma= as_scalar((sum(Sigma2,0)*(betayh - alpha*mu))/accu(Sigma2));
		 mat ww;
		 ww=ones<mat>(1,p2);
		 gamma=as_scalar((ww * (betayh*(n2-1) - alpha*Sigma2*mu)))/accu(Sigma2);
		 
      }   
      
      //ztmp = as_scalar(( (betayh - w2*gamma).t())*Sigma2*mu);
	  
	  ztmp = as_scalar(mu.t()*betayh*(n2-1)-gamma*mu.t()*sum(Sigma2,1));
      mat trtmp2 = Sigma2 * Sigma;
      double tr2 = trace(trtmp2);
      if (constralpha == 1){                                                     
        alpha = 0;
      }
      else {
        alpha = as_scalar(ztmp/((mu.t()*Sigma2*mu) + tr2));
      }
      alpha2 = alpha*alpha;
	  
	  
	  res_betaxh=as_scalar((n1-1)-2*(n1-1)*lambda*mu.t()*betaxh+lambda2*mu.t()*Sigma1*mu);
	    
      res_betayh = as_scalar((n2-1)-2*alpha*(n2-1)*mu.t()*betayh-2*gamma*(n2-1)*sum(betayh)+ alpha2*mu.t()*Sigma2*mu + gamma*gamma*accu(Sigma2)+2*alpha*gamma*mu.t()*sum(Sigma2,1));
	  
      
      //res_betaxh = betaxh - mu*lambda;
      //res_betayh = betayh - alpha*mu - w2*gamma; 
	 // cout<< mu << endl;
	   // cout<< lambda << endl;
	    //cout<< tr1 << endl;
      
      sigma2x = as_scalar((res_betaxh + lambda2*tr1)/n1);
	 // cout<< sigma2x << endl;
      sigma2y = as_scalar((res_betayh + alpha2*tr2)/n2);
	 // cout<< sigma2y << endl;
      sigma2beta = as_scalar((sum(mu % mu) + trace(Sigma))/p);
	 // cout<< sigma2beta << endl;
      
      // Reduction-step
      sigma2beta = lambda2 * sigma2beta;
      alpha = alpha / lambda;
      alpha2 = alpha*alpha;
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
	
	  if (p * sigma2beta < 1e-04){
    perror("The estimate of gene expression heritability explained by cis-SNPs is smaller than 0.01%");
  }
    
    
    return List::create(Rcpp::Named("alpha") = alpha,
                        Rcpp::Named("gamma") = gamma,
                        Rcpp::Named("sigmaX") = sigma2x,
                        Rcpp::Named("sigmaY") = sigma2y,
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
