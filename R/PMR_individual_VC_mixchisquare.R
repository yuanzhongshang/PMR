#' @title The main function for PMR model with individual-level data under polygenic pleiotropy effect assumption (variance component model),with mixed chisquare approximation for the test #' statistics
#' @description PMR_individual_VC_mixchisquare applies an EM argorithm to obtain the estimates of the variance component parameters and conduct a score test for causal effect test, with 
#' mixed chisquare approximation for the test statistics often when the sample size is small
#' @param yin standarized exposure vector (e.g. gene expression in TWAS)
#' @param zin standardized complex trait vector.
#' @param x1in standardized cis-genotype matrix in eQTL data.
#' @param x2in standardized cis-genotype matrix in GWAS data.
#' @param max_iterin The maximum iteration, which can be determined by users.
#' @return  the p value for the causal effect test
PMR_individual_VC_mixchisquare<-function(yin, zin, x1in, x2in, max_iterin=max_iterin){
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
inter=e2$values+sigma/sigmaa
intermatrix=e2$vectors%*%diag(1/inter)%*%t(e2$vectors)


HinverseAtilde<-(Atilde-Atilde%*%intermatrix%*%SigSig)/sigma
W=HinverseAtilde%*%M%*%t(A)
#observe###
q1=t(Y2)%*%W%*%(Y1)
svdW<-svd(W)
U=svdW$u
V=svdW$v
lamda<-diag(svdW$d)

tAtildeU<-t(Atilde)%*%U
sigma1<-t(tAtildeU)%*%tAtildeU
eigensigma1<-eigen(sigma1)

tAV<-t(A)%*%V
sigma2<-t(tAV)%*%tAV
eigensigma2<-eigen(sigma2)
sigma1values<-sigmaa*eigensigma1$values+sigma
sigma2values<-sigmag*eigensigma2$values+sigmae

sqrtsigma1<-eigensigma1$vectors%*%diag(sqrt(sigma1values))%*%t(eigensigma1$vectors)
sqrtsigma2<-eigensigma2$vectors%*%diag(sqrt(sigma2values))%*%t(eigensigma2$vectors)

Q=sqrtsigma1%*%lamda%*%sqrtsigma2
    s=sd(c(Q))
    e=try(svd(Q/s))
    lambda=c(0.5*e$d,-0.5*e$d)*s
   r1=davies(abs(q1), lambda, h = rep(1, length(lambda)), delta = rep(0, length(lambda)), sigma = 0, lim = 10000, acc = 0.0001)
    pvalue=2*r1$Qq
 return(pvalue)
}
	