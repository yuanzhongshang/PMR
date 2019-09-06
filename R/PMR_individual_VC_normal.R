#' @title The main function for PMR model with individual-level data under polygenic pleiotropy effect assumption (variance component model),with normal approximation for the test 
#' statistics
#' @description PMR_individual_VC_normal applies an EM argorithm to obtain the estimates of the variance component parameters and conduct a score test for causal effect test, with 
#' normal approximation for the test statistics often when the sample size is large
#' @param yin standarized exposure vector (e.g. gene expression in TWAS)
#' @param zin standardized complex trait vector.
#' @param x1in standardized cis-genotype matrix in eQTL data.
#' @param x2in standardized cis-genotype matrix in GWAS data.
#' @param max_iterin The maximum iteration, which can be determined by users.
#' @return the p value for the causal effect test 
PMR_individual_VC_normal<-function(yin, zin, x1in, x2in, max_iterin=max_iterin){

A<-x1in
Atilde<-x2in
Y2<-zin
Y1<-yin

n1<-nrow(A)
n2<-nrow(Atilde)
p<-ncol(A)
w2 = matrix(rep(1,n2),ncol=1)
w1 = matrix(rep(1,n1),ncol=1)
if(n1<=p){
fm0 = Estimate_VC_1(Y1,w1,A, max_iterin)
sigmag =as.vector(fm0$sigmaVC)
sigmae =as.vector(fm0$sigmaerror)
}else{
fm0 = Estimate_VC_2(Y1,w1,A,max_iterin)
sigmag =as.vector(fm0$sigmaVC)
sigmae =as.vector(fm0$sigmaerror)
}
if(n2<=p){
Yfm0 = Estimate_VC_1(Y2,w2,Atilde,max_iterin)
sigmaa =as.vector(Yfm0$sigmaVC)
sigma =as.vector(Yfm0$sigmaerror)
}else{
Yfm0 = Estimate_VC_2(Y2,w2,Atilde,max_iterin)
sigmaa =as.vector(Yfm0$sigmaVC)
sigma =as.vector(Yfm0$sigmaerror)
}
Sig<-t(A)%*%A
e1<-eigen(Sig)
SigSig<-t(Atilde)%*%Atilde
e2<-eigen(SigSig)
b=e1$values+sigmae/sigmag
M=e1$vectors%*%diag(1/b)%*%t(e1$vectors)
inter<-sigmaa/(e2$values*sigmaa+sigma)
intermatrix=e2$vectors%*%diag(inter)%*%t(e2$vectors)
tyAtilde<-t(Y2)%*%Atilde;
q1<-(tyAtilde-tyAtilde%*%intermatrix%*%SigSig)%*%M%*%t(A)%*%Y1/sigma;
observe<-as.vector(q1);
############################
tAtildeHinverseAtilde=(SigSig-SigSig%*%intermatrix%*%SigSig)/sigma;
Q=A%*%M%*%tAtildeHinverseAtilde%*%M%*%t(A)%*%(sigmag*A%*%t(A)+sigmae*diag(n1));
traceQ<-sum(diag(Q))
statistics<-observe/sqrt(traceQ)
pvalue<-2*(1-pnorm(abs(statistics)))
return(pvalue)
}
	