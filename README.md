# PMR

PMR(Probabilistic two sample mendelian randomization),is an R package for efficient statistical inference of two-sample MR analysis in transcriptomewide association studies (TWAS). It can account for the correlatded instruments and the horizontal pleiotropy, and can provide the
accurate estimates of both causal effect and horizontal pleiotropy effect as well as the two corresponding p values. PMR can be applied in single trait analysis as well as multiple correlated outcome traits analysis.

# Installation
It is easy to install the development version of PMR package using the 'devtools' package. The typical install time on a "normal" desktop computer is less than one minute.

```
# install.packages("devtools")
library(devtools)
install_github("yuanzhongshang/PMR")
```


# Usage
There are four main functions in PMR package, one is *PMR_individual* for single trait with individual level data, another is *PMR_summary_Egger* for single trait with summary statistics. For PMR_indvidual, two pleiotropy model assumptions have been designed, one is the Egger assumption, the other is polygenic assumption (variance component model). 
Note that the current version of pleiotropy assumption for summary data is Egger, which has been implemented by function *PMR_summary_Egger*. The other two corresponding functions are *moPMR_individual* for multiple outcome traits with individual level data, and *moPMR_summary* for multiple outcomes traits with summary statistics.
You can find the instructions by '?PMR_individual' , '?PMR_summary_Egger', '?moPMR_individual' and '?moPMR_summary'.
```
library(PMR)

?PMR_individual

?PMR_summary_Egger

?moPMR_individual

?moPMR_summary
```

# Example

One simple example to use the package can be found at https://github.com/yuanzhongshang/PMR/tree/master/example

# Results reproduced 

All results from all methods used in the PMR paper can be reproduced at https://github.com/yuanzhongshang/PMRreproduce for single trait analysis, and at https://github.com/yuanzhongshang/moPMRreproduce for multiple outcome traits analysis.

# Development
This R package is developed by Lu Liu, Zhongshang Yuan and Xiang Zhou.

