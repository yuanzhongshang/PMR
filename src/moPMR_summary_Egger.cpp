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
  //loglik:最大似???  datum::inf:???, infinity
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
    
    
    
    void loglike_twas(const arma::mat& R, const arma::vec& res_betaxh, const arma::vec& res_betayh, const arma::vec& mu, const double& sigma2z, const double& sigma2x,
    const arma::vec& eVal, const int& n1, const int& n2, const int& p, double& loglik){
    double Eab;//, loglik;
    Eab = as_scalar(res_betaxh)/2/sigma2x + as_scalar(res_betayh)/2 + sum(mu % mu)/2/sigma2z;
    //cout<<res_betayh<<endl;
    loglik = -log(sigma2z)*p*0.5 - log(sigma2x)*n1*0.5 - accu(log(eVal))*n2*0.5 - sum(log(R.diag()))- Eab;
    //cout<<eVal(find(eVal > 0.00001))<<endl;
    }
    
    
    
    // [[Rcpp::export]]
    SEXP moPMR_summary_Egger(SEXP betaxin, SEXP betayin, SEXP Sigma1sin, SEXP Sigma2sin, SEXP Sigmayin, SEXP samplen1, SEXP samplen2, SEXP gammain, SEXP alphain, SEXP max_iterin, SEXP epsin){// *
    try{
    const int constrgamma = Rcpp::as<int>(gammain);
    const int constralpha = Rcpp::as<int>(alphain);
    const int n1 = Rcpp::as<int>(samplen1);  
    const int n2 = Rcpp::as<int>(samplen2); 
    const int maxIter = Rcpp::as<int>(max_iterin);
    const double epsStopLogLik = Rcpp::as<double>(epsin);
    const int pxem_indicator = 1;
    
    const arma::vec betax = as<arma::vec>(betaxin);
    const arma::mat betay = as<arma::mat>(betayin);
    const arma::mat Sigma1s = as<arma::mat>(Sigma1sin);
    const arma::mat Sigma2s = as<arma::mat>(Sigma2sin);
    const arma::mat Sigmays = as<arma::mat>(Sigmayin);
    //const arma::vec constralpha = as<arma::vec>(alphain);
    // *convert data format
    int p1 = Sigma1s.n_cols, p2 = Sigma2s.n_cols, k = betay.n_rows;
    mat w2;
    
    if (p1 != p2){
    perror("The dimensions of Zx and Zy are not matched");
    }
    int p = p1;
    if (constrgamma == 1){
    w2 = zeros(p2,1);
    }
    else {
    w2 = ones(p2,1);
    }
    //cout << w2 << endl;
    vec betaxh=betax;
    mat betayh=betay;
    //cout << betayh << endl;
    mat Sigma1 = Sigma1s* (n1-1);
    mat Sigma2 = Sigma2s* (n2-1);
    mat Sigmay = Sigmays* (n2-1);
    
    //double w2tw2 = as_scalar(w2*w2.t());
    // initialization using lmm_pxem
    double sigma2x, sigma2z, loglik0;
    
    int iter0;
    mat Sigb = zeros<mat>(p,p);
    vec mub  = zeros<vec>(p);
    
    lmm_pxem_ptr2_summary(betaxh, Sigma1, n1, maxIter, sigma2x, sigma2z, loglik0, iter0, Sigb, mub);
    
    //declaration of variables used within loop
    mat Sigma_inv(p,p), R(p,p),invR(p,p), Sigma(p,p), res_betayh1(k,k);
    vec mutmp(p), mu(p), res_betaxh(p), res_betayh(p), gamma(k);
    mat I = eye<mat>(p,p);
    //cout << mean(w2) << endl;
    //initialization of parameters
    if (constrgamma == 1){
    gamma = zeros(k,1);
    }
    else {
    gamma =  betayh * w2 * ((n2-1)/ accu(Sigma2));
    }
    //cout << gamma << endl;
    // vec gamma = solve(w2tw2, w2ty);
    double lambda = 1, lambda2 = lambda*lambda;
    mat  sigmay = Sigmays, S(k,k),sigmay_inv(k,k);
    vec alpha = zeros<vec>(k,1);
    //cout << sigmay << endl;
    // vec eVal;
    // mat eVec;
    // eig_sym(eVal, eVec, sigmay);
    //cout << eVec << endl;
    // S = diagmat(1/eVal);
    //cout << S << endl;
    // sigmay_inv = eVec*S*eVec.t();
    //cout << sigmay_inv << endl;
    
    //cout << ZytZy_inv << endl;
    vec loglik(maxIter), loglik2(maxIter);
    
    // initialization of likelihood
    loglik(0) = NAN;
    int Iteration = 1;
    for (int iter = 2; iter <= maxIter; iter ++ ) {
    // E-step
    vec eVal;
    mat eVec;
    eig_sym(eVal, eVec, sigmay);
    //cout << eVal << endl;
    S = diagmat(1.0/eVal);
    //cout << S << endl;
    sigmay_inv = eVec*S*eVec.t();
    Sigma_inv = Sigma1/sigma2x + Sigma2*as_scalar(alpha.t()*sigmay_inv*alpha) + I /sigma2z;
    //cout << Sigma_inv << endl;
    // vec eVal1;
    // mat eVec1;
    // eig_sym(eVal1, eVec1, Sigma_inv);
    // //cout << eVal1 << endl;
    // R = diagmat(1.0/eVal1);
    // Sigma = eVec1*R*eVec1.t();
    R = chol(Sigma_inv);
    invR = inv(R);
    Sigma = invR*invR.t();
    //cout << Sigma << endl;
    
    // mat L= zeros<mat>(p,n2);
    // L.each_row() += alpha.t()*sigmay_inv*(Y - gamma*w2);
    // cout << L << endl;
    // muty = sum(Zy.t() % L,1);
    //cout << muty << endl;
    mutmp = (n1-1)*betaxh/sigma2x + ((n2-1)*betayh - gamma*sum(Sigma2,0)).t()*(sigmay_inv*alpha);
    mu = Sigma*mutmp;
    //cout << mu(1) << endl;
    //evaluate incomplete log-likelihood
    res_betaxh=as_scalar((n1-1)-2*(n1-1)*lambda*mu.t()*betaxh+lambda2*mu.t()*Sigma1*mu);
    //cout << res_betaxh << endl;
    //res_betayh = as_scalar((n2-1)-2*alpha*(n2-1)*mu.t()*betayh-2*gamma*(n2-1)*sum(betayh)+ alpha2*mu.t()*Sigma2*mu + gamma*gamma*accu(Sigma2)+2*alpha*gamma*mu.t()*sum(Sigma2,1));
    res_betayh = as_scalar(trace(sigmay_inv*Sigmays*(n2-1))-2*(n2-1)*trace(sigmay_inv*betayh*(mu*alpha.t()+w2*gamma.t()))+as_scalar(alpha.t()*sigmay_inv*alpha*mu.t()*Sigma2*mu)+as_scalar(gamma.t()*sigmay_inv*alpha*mu.t()*Sigma2*w2)+as_scalar(alpha.t()*sigmay_inv*gamma*w2.t()*Sigma2*mu)+as_scalar(gamma.t()*sigmay_inv*gamma)*accu(Sigma2));
    //mat s= sigmay_inv*Sigmays;
    //mat ss= sigmay_inv*betayh*(mu*alpha.t()+w2*gamma.t());
    //res_betayh = as_scalar(sum(s.diag())*(n2-1))-2*(n2-1)*sum(ss.diag())+as_scalar(alpha.t()*sigmay_inv*alpha*mu.t()*Sigma2*mu)+as_scalar(gamma.t()*sigmay_inv*alpha*mu.t()*Sigma2*w2)+as_scalar(alpha.t()*sigmay_inv*gamma*w2.t()*Sigma2*mu)+as_scalar(gamma.t()*sigmay_inv*gamma)*accu(Sigma2);
    //res_betayh = trace(sigmay_inv*(n2-1)*(Sigmays-betayh*mu*alpha.t()-betayh*w2*gamma.t()-alpha*mu.t()*betayh.t()-gamma*w2.t()*betayh.t())+as_scalar(mu.t()*Sigma2*mu)*alpha*alpha.t()+as_scalar(mu.t()*Sigma2*w2)*alpha*gamma.t()+as_scalar(w2.t()*Sigma2*mu)*gamma*alpha.t()+gamma*gamma.t()*accu(Sigma2));
    //cout << res_betayh << endl;
    //double res_betayhsum=0;
    
    // for(l=0;l<n2;l++)
    // res_ysum = res_ysum + as_scalar((Y.col(l) - alpha*(Zymu.row(l)) -gamma*(w2.col(l))).t()*sigmay_inv*(Y.col(l) - alpha*(Zymu.row(l)) -gamma*(w2.col(l))));
    //mat s= res_betayh.t()*sigmay_inv*res_betayh;
    //res_betayhsum = sum(s.diag());
    //loglik(iter - 1) = loglike_twas(R, res_x, res_y, mu, sigma2z, sigma2x, sigma2y,n1,n2,p);
    loglike_twas(R, res_betaxh, res_betayh, mu, sigma2z, sigma2x, eVal, n1 , n2, p, loglik(iter - 1));
    
    if ( loglik(iter - 1) - loglik(iter - 2) < -1e-7){
    perror("The likelihood failed to increase!");
    }
    // M-step
    mat trtmp1 = Sigma1  * Sigma;
    double tr1 = trace(trtmp1);
    
    if (pxem_indicator == 1){
    lambda = as_scalar(( betaxh.t() *(n1-1)* mu )/((mu.t() * Sigma1 * mu) + tr1));
    lambda2 = lambda*lambda;
    }
    else {
    lambda = 1; lambda2 = 1;
    }
    
    //mux = solve(w1tw1, w1tx - lambda*Zxtw1.t()*mu);
    //cout << lambda << endl;
    
    //mat gamma = zeros<mat>(k,1);
    if (constrgamma == 1){
    // w2= ones(n2,1);
    gamma = zeros(k,1);
    }
    else {
    gamma = (betayh*(n2-1)-alpha*mu.t()*Sigma2)*w2/accu(Sigma2);
    }
    //cout << gamma << endl;
    //gamma= solve(w2tw2, w2ty - alpha*(Zytw2.t()*mu));
    // vec ytmp = zeros<vec>(k);
    // for(l=0;l<n2;l++)
    // ytmp = as_scalar(Zymu.row(l))*(Y.col(l) - gamma*w2);
    mat trtmp2 = Sigma2  * Sigma;
    double tr2 = trace(trtmp2);
    
    
    //uvec q =  find(constralpha == 0);
    //cout << q <<endl;
    // int i;
    // for(i=0;i<k-accu(constralpha);i++){
    // 
    //   alpha(q(i)) = as_scalar(((Y.row(q(i)) - gamma.row(q(i))*w2) * Zymu))/(sum(Zymu % Zymu) + tr2);
    // 
    // }
    
    //alpha(q) = (((n2-1)*betayh.rows(q) - gamma.rows(q)*sum(Sigma2,0)) * mu)/((mu.t()*Sigma2*mu) + tr2);
    if (constralpha == 1){
      mat alpha = zeros<mat>(k,1);
    }
    else {
      alpha = (((n2-1)*betayh - gamma*sum(Sigma2,0)) * mu)/as_scalar((mu.t()*Sigma2*mu) + tr2);
    }
    //cout << alpha <<endl;
    // alpha.each_row(find(constralpha == 0)) = 0;
    // alpha.each_row(find(constralpha)) = (Y -gamma*w2)(find(constralpha))* Zymu/(sum(Zymu % Zymu) + tr2);
    
    mat alpha2(p,p);
    alpha2 = alpha*(alpha.t());
    //cout << alpha2 <<endl;
    res_betaxh=as_scalar((n1-1)-2*(n1-1)*lambda*mu.t()*betaxh+lambda2*mu.t()*Sigma1*mu);
    res_betayh1 = (n2-1)*(Sigmays-betayh*mu*alpha.t()-betayh*w2*gamma.t()-alpha*mu.t()*betayh.t()-gamma*w2.t()*betayh.t())+as_scalar(mu.t()*Sigma2*mu)*alpha*alpha.t()+as_scalar(mu.t()*Sigma2*w2)*alpha*gamma.t()+as_scalar(w2.t()*Sigma2*mu)*gamma*alpha.t()+gamma*gamma.t()*accu(Sigma2);
    //cout << res_betayh1 <<endl;
    sigma2x = as_scalar((res_betaxh + lambda2*tr1)/n1);
    sigmay = (res_betayh1 + alpha2 * tr2)/n2;
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