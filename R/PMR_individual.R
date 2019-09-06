#' @title The main function for probabilistic two-sample MR model with individual-level data
#' @description PMR_individual applies a likelihood-based approach, accounting for the correlated instruments and horizontal pleiotropy effect 
#' @param yin standarized exposure vector (e.g. gene expression in TWAS)
#' @param zin standardized complex trait vector.
#' @param x1in standardized cis-genotype matrix in eQTL data.
#' @param x2in standardized cis-genotype matrix in GWAS data.
#' @param method indicate the methods with Egger pleiotropy effect assumption (method="PMR_individual_Egger"), polygenic pleiotropy effect assumption with mixed chisquare approximation for the test statistics(method="PMR_individual_VC_mixchisquare") when the sample size is small or polygenic pleiotropy effect assumption with normal approximation for the test statistics(method="PMR_individual_VC_normal") when the sample size is large.
#' @param max_iterin The maximum iteration, which can be determined by users.
#' @param epsin The convergence tolerance of the absolute value of the difference  between the nth and (n+1)th log likelihood, which can be determined by users.
#' @param Heritability_geneexpression_threshold The threshold for the estimate of gene expression heritability explained by cis-SNPs, which can be determined by users. The causal effect pvalue will be assigned to be NA automatically if the the estimate of gene expression heritability is under this threshold.
#' @return A list of estimated parameters including the p values for both the causal effect test and pleiotropy effect test. Note that the current output under polygenic pleiotropy effect assumption is only causal effect p values. 
#' \item{causal_effect}{The estimates of causal effect}
#' \item{causal_pvalue}{The p value for the causal effect}
#' \item{pleiotropy_effect}{The estimates of pleiotropy effect}
#' \item{pleiotropy_pvalue}{The p value for the pleiotropy effect}
#' \item{sigma_cisSNP}{The variance estimate for the cis-SNP effect sizes}
#' \item{sigma_error_1}{The variance estimate of the error in eQTL data model}
#' \item{sigma_error_2}{The variance estimate of the error in GWAS data model}
PMR_individual<-function(yin, zin, x1in, x2in, method="PMR_individual_Egger", max_iterin =1000,epsin=1e-5, Heritability_geneexpression_threshold=1e-04){

eQTLdata<-na.omit(cbind(yin,x1in))
yin_NA<-as.vector(eQTLdata[,1])
x1in_NA<-eQTLdata[,-1]

GWASdata<-na.omit(cbind(zin,x2in))
zin_NA<-as.vector(GWASdata[,1])
x2in_NA<-GWASdata[,-1]

x1in_NA_sd<-as.vector(apply(x1in_NA,2,sd))

x2in_NA_sd<-as.vector(apply(x2in_NA,2,sd))

if(sd(yin_NA)!=0 & sd(zin_NA)!=0 & sum(x1in_NA_sd==0)==0 & sum(x2in_NA_sd==0)==0){

y<-as.vector(scale(yin_NA))
z<-as.vector(scale(zin_NA))
x1<-scale(x1in_NA)
x2<-scale(x2in_NA)

#################################################################################
if(method=="PMR_individual_Egger"){
fmH1 = PMR_individual_Egger(y, z, x1, x2,gammain=0,alphain = 0, max_iterin=max_iterin,epsin=epsin)
p=ncol(x1)
Heritability_estimate=p*fmH1$sigmabeta

##test whether the gene expression heritability explained by cis-SNPs is small enough

if(Heritability_estimate <= Heritability_geneexpression_threshold){
cat(paste0("## The estimate of gene expression heritability explained by cis-SNPs is smaller than the threshold ##") )
cat(paste0("## The pvalue of the causal effect test was assigned to be NA ##") )
fmH0gamma = PMR_individual_Egger(y, z, x1, x2, gammain=1, alphain = 0, max_iterin=max_iterin,epsin=epsin)
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
fmH0gamma = PMR_individual_Egger(y, z, x1, x2, gammain=1, alphain = 0, max_iterin=max_iterin,epsin=epsin)
fmH0alpha = PMR_individual_Egger(y, z, x1, x2, gammain=0, alphain = 1, max_iterin=max_iterin,epsin=epsin)
loglikH1=max(fmH1$loglik,na.rm=T)
loglikH0gamma=max(fmH0gamma$loglik,na.rm=T)
loglikH0alpha=max(fmH0alpha$loglik,na.rm=T)
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
if(method=="PMR_individual_VC_mixchisquare"){
pvalue<-PMR_individual_VC_mixchisquare(y, z, x1, x2, max_iterin=max_iterin)
result=list()
result$causal_pvalue<-pvalue
return(result)
}
if(method=="PMR_individual_VC_normal"){
pvalue<-PMR_individual_VC_normal(y, z, x1, x2, max_iterin=max_iterin)
result=list()
result$causal_pvalue<-pvalue
return(result)
}
######################
}
}




