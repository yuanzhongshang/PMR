#' @title The main function for probabilistic two-sample MR model for multiple outcome traits with individual-level data
#' @description moPMR_individual applies a likelihood-based approach, accounting for the correlated instruments and horizontal pleiotropy effect 
#' @param xin standarized exposure vector (e.g. gene expression in TWAS)
#' @param Yin complex trait matrix, each column is the standardized specific trait.
#' @param Zxin standardized cis-genotype matrix in eQTL data.
#' @param Zyin standardized cis-genotype matrix in GWAS data.
#' @param max_iterin The maximum iteration, which can be determined by users.
#' @param epsin The convergence tolerance of the absolute value of the difference  between the nth and (n+1)th log likelihood, which can be determined by users.
#' @param Heritability_geneexpression_threshold The threshold for the estimate of gene expression heritability explained by cis-SNPs, which can be determined by users. The causal effect pvalue will be assigned to be NA automatically if the the estimate of gene expression heritability is under this threshold.
#' @return A list of estimated parameters including the p values for both the causal effect test and pleiotropy effect test. 
#' \item{causal_effect}{The estimates of causal effect}
#' \item{causal_pvalue}{The p value for the causal effect}
#' \item{pleiotropy_effect}{The estimates of pleiotropy effect}
#' \item{pleiotropy_pvalue}{The p value for the pleiotropy effect}
#' \item{sigma_cisSNP}{The variance estimate for the cis-SNP effect sizes}
#' \item{sigma_error_1}{The variance estimate of the error in eQTL data model}
#' \item{sigma_error_2}{The variance estimate of the error in GWAS data model}
moPMR_individual<-function(xin, Yin, Zxin, Zyin, max_iterin =1000,epsin=1e-5, Heritability_geneexpression_threshold=1e-04){
  
  eQTLdata<-na.omit(cbind(xin,Zxin))
  xin_NA<-as.vector(eQTLdata[,1])
  Zxin_NA<-eQTLdata[,-1]
  
  k=dim(Yin)[2]
  GWASdata<-na.omit(cbind(Yin,Zyin))
  Yin_NA<-as.matrix(GWASdata[,1:k])
  Zyin_NA<-GWASdata[,-c(1:k)]
  
  Zxin_NA_sd<-as.vector(apply(Zxin_NA,2,sd))
  
  Zyin_NA_sd<-as.vector(apply(Zyin_NA,2,sd))
  
  if(sd(Zxin_NA)!=0 & sd(Zyin_NA)!=0 & sum(Zxin_NA_sd==0)==0 & sum(Zyin_NA_sd==0)==0){
    
    x<-as.vector(scale(xin_NA))
    Y<-as.matrix(t(scale(Yin_NA)))
    Zx<-scale(Zxin_NA)
    Zy<-scale(Zyin_NA)
    
    #################################################################################
    
    fmH1 = moPMR_individual_Egger(x, Y, Zx, Zy,gammain=0,alphain = 0, max_iterin=max_iterin,epsin=epsin)
    p=ncol(Zx)
    Heritability_estimate=p*fmH1$sigmaZ
    
    ##test whether the gene expression heritability explained by cis-SNPs is small enough
    
    if(Heritability_estimate <= Heritability_geneexpression_threshold){
      cat(paste0("## The estimate of gene expression heritability explained by cis-SNPs is smaller than the threshold ##") )
      cat(paste0("## The pvalue of the causal effect test was assigned to be NA ##") )
      fmH0gamma = moPMR_individual_Egger(x, Y, Zx, Zy, gammain=1, alphain = 0, max_iterin=max_iterin,epsin=epsin)
      loglikH1=max(fmH1$loglik,na.rm=T)
      loglikH0gamma=max(fmH0gamma$loglik,na.rm=T)
      stat_gamma = 2 * (loglikH1 - loglikH0gamma)
      pvalue_gamma = pchisq(stat_gamma,k,lower.tail=F)
      result=list()
      result$causal_effect=fmH1$alpha
      result$causal_pvalue=NA
      result$pleiotropy_effect=fmH1$gamma
      result$pleiotropy_pvalue=pvalue_gamma
      result$sigma_cisSNP=fmH1$sigmaZ
      result$sigma_error_1=fmH1$sigmaX
      result$sigma_error_2=fmH1$sigmaY
      return(result)
    }else{
      fmH0gamma = moPMR_individual_Egger(x, Y, Zx, Zy, gammain=1, alphain = 0, max_iterin=max_iterin,epsin=epsin)
      fmH0alpha = moPMR_individual_Egger(x, Y, Zx, Zy, gammain=0, alphain = 1, max_iterin=max_iterin,epsin=epsin)
      loglikH1=max(fmH1$loglik,na.rm=T)
      loglikH0gamma=max(fmH0gamma$loglik,na.rm=T)
      loglikH0alpha=max(fmH0alpha$loglik,na.rm=T)
      stat_alpha = 2 * (loglikH1 - loglikH0alpha)
      pvalue_alpha = pchisq(stat_alpha,k,lower.tail=F)
      stat_gamma = 2 * (loglikH1 - loglikH0gamma)
      pvalue_gamma = pchisq(stat_gamma,k,lower.tail=F)
      result=list()
      result$causal_effect=fmH1$alpha
      result$causal_pvalue=pvalue_alpha
      result$pleiotropy_effect=fmH1$gamma
      result$pleiotropy_pvalue=pvalue_gamma
      result$sigma_cisSNP=fmH1$sigmaZ
      result$sigma_error_1=fmH1$sigmaX
      result$sigma_error_2=fmH1$sigmaY
      return(result)
    }
    
    ######################
  }
}
