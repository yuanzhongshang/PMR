#' @importFrom PDSCE pdsoft
#' @title The main function for probabilistic two-sample MR model with summary data under Egger pleiotropy effect assumption
#' @description PMR_summary_Egger applies a likelihood-based approach, accounting for the correlated instruments and horizontal pleiotropy effect 
#' @param Zscore_1 the Zscore vector of the cis-SNP effect size vector for one specific gene in eQTL data
#' @param Zscore_2 the Zscore vector of the cis-SNP effect size vector for one specific gene in GWAS data
#' @param Sigma1sin the LD matrix in eQTL data
#' @param Sigma2sin the LD matrix in GWAS data,both Sigma2sin and sigma1sin are often the same from the reference panel
#' @param samplen1 the sample size of eQTL data
#' @param samplen2 the sample size of GWAS data
#' @param lambda The shrinkage parameter to guarantee the sparsity and positive definiteness of the estimated GWAS LD matrix using
#' the method in PDSCE package,the default value is 0 to indicate almost no shrikage, other choice can be empirically 
#' chosen to be 0.05,0.1,0.15. 
#' @param max_iterin The maximum iteration, which can be determined by users.
#' @param epsin The convergence tolerance of the absolute value of the difference  between the nth and (n+1)th log likelihood, which can be determined by users.
#' @param Heritability_geneexpression_threshold The threshold for the estimate of gene expression heritability explained by cis-SNPs, which can be determined by users. The causal effect pvalue will be assigned to be NA automatically if the the estimate of gene expression heritability is under this threshold.
#' @return A list of estimated parameters including the p values for both the causal effect test and pleiotropy effect test
#' \item{causal_effect}{The estimates of causal effect}
#' \item{causal_pvalue}{The p value for the causal effect}
#' \item{pleiotropy_effect}{The estimates of pleiotropy effect}
#' \item{pleiotropy_pvalue}{The p value for the pleiotropy effect}
#' \item{sigma_cisSNP}{The variance estimate for the cis-SNP effect sizes}
#' \item{sigma_error_1}{The variance estimate of the error in eQTL data model}
#' \item{sigma_error_2}{The variance estimate of the error in GWAS data model}

PMR_summary_Egger<-function(Zscore_1, Zscore_2, Sigma1sin, Sigma2sin, samplen1, samplen2, lambda=0, max_iterin =1000,epsin=1e-5, Heritability_geneexpression_threshold=1e-04){
betaxin<-Zscore_1/sqrt(samplen1-1)
betayin<-Zscore_2/sqrt(samplen2-1)
Sigma2sin<-pdsoft(Sigma2sin,lam=lambda)$theta 
fmH1=PMR_summary_Egger_CPP(betaxin,betayin,Sigma1sin,Sigma2sin,samplen1,samplen2,gammain=0,alphain=0,max_iterin =max_iterin, epsin=epsin)
p=length(betaxin)
Heritability_estimate=p*fmH1$sigmabeta
##test whether the gene expression heritability explained by cis-SNPs is small enough
if(Heritability_estimate <= Heritability_geneexpression_threshold){
cat(paste0("## The estimate of gene expression heritability explained by cis-SNPs is smaller than the threshold ##") )
cat(paste0("## The pvalue of the causal effect test was assigned to be NA ##") )
fmH0gamma=PMR_summary_Egger_CPP(betaxin, betayin,Sigma1sin, Sigma2sin, samplen1, samplen2, gammain=1,alphain=0, max_iterin =max_iterin, epsin=epsin)
loglikH1=max(fmH1$loglik,na.rm=T)
loglikH0gamma=max(fmH0gamma$loglik,na.rm=T)
stat_gamma = 2 * (loglikH1 - loglikH0gamma)
pvalue_gamma = pchisq(stat_gamma,1,lower.tail=F)
result=list()
result$causal_effect=fmH1$alpha
result$causal_pvalue=NA
result$pleiotropy_effect=fmH1$gamma
result$pleiotropy_pvalue=pvalue_gamma
result$sigma_cisSNP=fmH1$sigmabeta
result$sigma_error_1=fmH1$sigmaX
result$sigma_error_2=fmH1$sigmaY
return(result)
}else{
fmH0alpha=PMR_summary_Egger_CPP(betaxin, betayin, Sigma1sin, Sigma2sin, samplen1, samplen2, gammain=0, alphain=1,max_iterin =max_iterin, epsin=epsin)
fmH0gamma=PMR_summary_Egger_CPP(betaxin, betayin,Sigma1sin, Sigma2sin, samplen1, samplen2, gammain=1,alphain=0, max_iterin =max_iterin, epsin=epsin)
loglikH1=max(fmH1$loglik,na.rm=T)
loglikH0alpha=max(fmH0alpha$loglik,na.rm=T)
loglikH0gamma=max(fmH0gamma$loglik,na.rm=T)
stat_alpha = 2 * (loglikH1 - loglikH0alpha)
pvalue_alpha = pchisq(stat_alpha,1,lower.tail=F)
stat_gamma = 2 * (loglikH1 - loglikH0gamma)
pvalue_gamma = pchisq(stat_gamma,1,lower.tail=F)
result=list()
result$causal_effect=fmH1$alpha
result$causal_pvalue=pvalue_alpha
result$pleiotropy_effect=fmH1$gamma
result$pleiotropy_pvalue=pvalue_gamma
result$sigma_cisSNP=fmH1$sigmabeta
result$sigma_error_1=fmH1$sigmaX
result$sigma_error_2=fmH1$sigmaY
return(result)
}
}
