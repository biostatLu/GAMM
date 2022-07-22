# GAMM：GRS-informed gene-based association mixed model
# Background
Although several integrative methods have been proposed, how to incorporate trans-ethnic genetic risk scores (GRS) more efficiently remains less understood. The mixed effects score test (MiST) was developed for SNP-set association studies by modeling impacts of genetic variants based on functional annotations while allowing for variant-specific influence. However, MiST conducts score-based test, thus cannot quantify the relative contribution of trans-ethnic GRS and genetic loci to phenotypic variation. Instead, such associations might be driven by cis-SNPs alone, which also partly explains why many identified genes in previous GWASs were generally located near or within clusters of associated GWAS loci. 

To integrate trans-ethnic GRS into GWAS, we here proposed a novel statistical method, called GAMM (GRS-informed gene-based association mixed model), by modeling the effects of cis-SNPs of a given gene as a function of trans-ethnic GRS. As would be shown, GAMM formulates the integrative method (e.g., MiST) and individual-level methods (e.g., SKAT) within a unified framework through mixed models and includes many prior methods/tests as special cases. Under the context of polygenic genetic architecture that has been widely confirmed for many complex traits, we estimate α in terms of their marginal estimates while accounting for correlation among SNPs with a linkage disequilibrium (LD) matrix.The marginal estimates of αq are directly obtained from summary statistics of the trait in the base population, and LD is calculated from a reference panel that matches the base population. To efficiently estimate unknown parameters in GAMM, we developed a parameter expansion expectation maximum (PX-EM) algorithm that can be scalable to large-scale biobank data. We defined two useful quantities to measure relative genetic contributions of trans-ethnic GRS and its cis-SNPs to phenotypic variance. In addition, a fast and powerful likelihood ratio test method was proposed in GAMM to jointly test the total effects of cis-SNPs and trans-ethnic GRS on the phenotype. With extensive simulation studies, we demonstrated that GAMM can correctly maintain the type I error rate control when jointly testing the total genetic effects and was often more powerful to identify true association signals across various scenarios compared to existing methods.

GAMM is implemented in R statistical environment.
# Example
For joint effect test using LRT in GAMM
```ruby
library(data.table)
library(glmnet)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(doParallel)
sourceCpp("lmm_pxem.cpp")
sourceCpp("lmm_pxem2_ZPcisSNP_variance_NotationPXEM_EM.cpp")
source("LRTSim_lmm_PXEM_Rcpp.R")
source("estimate_beta.R")
weight=estimate_beta(BETA,SE,R)
g = as.matrix(apply(G%*%weight,2,scale)[,1])
fit = lmm_pxem2_ZPcisSNP(y, X=cbind(1, X, g), G=G, PXEM=TRUE, maxIter=1000)
simLike <- LRTSimZP(Z = X, E = g, G = G, nsim=1e5, parallel=c("multicore"), ncpus = 1L) ## exact LRT

fitx=lm(y~X)
obsLike = c((fit$loglik - logLik(fitx))*2)
plrt = mean(simLike >= obsLike)
fitmix = pmixLRT(simLike, c(1e3,1e4,1e5)) ## approximate LRT

//' @param weight  true effects for SNPs available for the same trait in another base population
//' @param y  response variable for GWAS data
//' @param X  covariates for GWAS data, not include the vector of ones
//' @param g  GRS = G*weight is the trans-enthic GRS information
//' @param G  genotype (cis-SNPs) matrix for GWAS

```

# Cite
Haojie Lu<sup>$</sup>, Shuiping Huang and Ping Zeng<sup>#</sup> (2021). A flexible and powerful gene-based association analysis framework for complex phenotypes in minority populations by integrating trans-ethnic genetic risk scores.

# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn.

# Update
2022-07-22 GAMM version 1.0
