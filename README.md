# GAMM：trans-ethnic GRS-informed gene-based association mixed model
# Background
Conducting association studies for traits in underrepresented populations even with much smaller sample sizes of participants than those in the EUR population is urgent. As GWAS sample sizes increase, the number of loci detected with statistical significance would also increase; therefore, the first solution is to recruit sufficient individuals from understudied populations, which could certainly improve power and discover more trait-associated loci. However, this is not feasible immediately because it is time-consuming and expensive. Alternatively, efficient statistical methods are developed to utilize existing data by incorporating widespread trans-ethnic genetic overlap shared by the same trait across diverse continental populations. Indeed, trans-ethnic genomics studies have been increasingly performed, and detected a great deal of new novel loci related to traits in non-EUR populations. An important finding of those multi-ancestry GWASs is that a large number of association signals identified in the EUR population are reproducible in understudied populations in the sense that significant loci exhibit high consistence in effect magnitude and association direction. However, how to leverage such shared information more efficiently in association analysis is less studied for traits in underrepresented populations.

In the present work, we refer to the primary population (e.g., EUR and/or EAS) as the base population, and the underrepresented population (e.g., AFR) as the target population. Our goal is to develop efficient association analysis strategies for complex traits in the target population with the hope of discovering more associated loci by leveraging genetic similarity between the target and base populations. To achieve this, we propose a novel statistical method called GAMM (trans-ethnic GRS-informed gene-based association mixed model), by modeling the effects of cis-SNPs of a given gene as a function of trans-ethnic genetic risk score (GRS). As would be demonstrated, GAMM formulates the integrative method of trans-ethnic GRS and individual-level gene-based methods (e.g., SKAT) within a unified framework through hierarchical mixed models and includes these two methods as special cases. Under the context of polygenic genetic architecture that has been widely confirmed for many complex traits, we estimate α in terms of their marginal estimates while accounting for correlation among SNPs with a linkage disequilibrium (LD) matrix.The marginal estimates of αq are directly obtained from summary statistics of the trait in the base population, and LD is calculated from a reference panel that matches the base population. To efficiently estimate unknown parameters in GAMM, we apply the parameter expansion expectation maximum (PX-EM) algorithm. Furthermore, a powerful likelihood ratio test (LRT) method is proposed to jointly test the total effects of SNPs and trans-ethnic GRS on the trait. We also define two useful quantities to quantity relative genetic contributions of the indirect trans-ethnic GRS effect and the direct SNP effect to the phenotypic variance. With extensive simulation studies, we demonstrate that GAMM could correctly maintain the type I error rate control when jointly examining the total genetic effects and is often more powerful to identify true association signals across various scenarios compared to competitive approaches.

GAMM is implemented in R statistical environment.
# Example
For the parameter estimation in GAMM
```ruby
library(glmnet)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("lmm_pxem2_ZPcisSNP_variance_NotationPXEM_EM.cpp")
source("LRTSim_lmm_PXEM_Rcpp.R")
source("estimate_beta.R")
BETA=read.table("BETA.txt",head=T) ## BETA is the effect size of SNPs in base population from summary statistics and SE is their standard error.
SE=read.table("SE.txt",head=T)
R=read.table("R.txt",head=T) ## R is the linkage disequilibrium (LD) matrix which can be calculated with genotypes of population-matched individuals from external reference panels such as the 1000 Genomes Project in the base population. We calculate R in a shrinkage fashion R = 0.95 * as.matrix(cor(G_base_population)) + diag(1-0.95, nrow=nSNPs, ncol=nSNPs)
G=read.table("G.txt",head=T)
y=read.table("y.txt",head=T)
X=read.table("X.txt",head=T)
BETA=as.matrix(BETA)
SE=as.matrix(SE)
R=as.matrix(R)
G=as.matrix(G)
y=as.matrix(y)
X=as.matrix(X)
weight=estimate_beta(BETA,SE,R)
g = as.matrix(apply(G%*%weight,2,scale)[,1])
fit = lmm_pxem2_ZPcisSNP(y, X=cbind(1, X, g), G=G, PXEM=TRUE, maxIter=1000)

v1 = var(g%*%fit$alpha[4,1])
v2 = var(G%*%fit$mub)
v3 = fit$sigma2e
v4 = var(cbind(1,X)%*%fit$alpha[1:3,1])
pve = (v1+v2)/(v1+v2+v3+v4)
pge = v1/(v1+v2)

//' @param y  response variable for GWAS data
//' @param X  covariates for GWAS data
//' @param g  GRS = G*weight is the trans-enthic GRS information
//' @param G  genotype (cis-SNPs) matrix for GWAS
//' @param maxIter  maximum iteration (default is 1000)

```
For joint effect test using LRT in GAMM
```ruby
library(glmnet)
library(Rcpp)
library(RcppArmadillo)
source("LRTSim_lmm_PXEM_Rcpp.R")
sourceCpp("lmm_pxem2_ZPcisSNP_variance_NotationPXEM_EM.cpp")
sourceCpp("LRTSim_lmm_PXEM_Rcpp.cpp")
source("estimate_beta.R")
BETA=read.table("BETA.txt",head=T) ## BETA is the effect size of SNPs in base population from summary statistics and SE is their standard error.
SE=read.table("SE.txt",head=T)
R=read.table("R.txt",head=T) ## R is the linkage disequilibrium (LD) matrix which can be calculated with genotypes of population-matched individuals from external reference panels such as the 1000 Genomes Project in the base population. We calculate R in a shrinkage fashion R = 0.95 * as.matrix(cor(G_base_population)) + diag(1-0.95, nrow=nSNPs, ncol=nSNPs)
G=read.table("G.txt",head=T)
y=read.table("y.txt",head=T)
X=read.table("X.txt",head=T)
BETA=as.matrix(BETA)
SE=as.matrix(SE)
R=as.matrix(R)
G=as.matrix(G)
y=as.matrix(y)
X=as.matrix(X)
weight=estimate_beta(BETA,SE,R)
g = as.matrix(apply(G%*%weight,2,scale)[,1])
fit = lmm_pxem2_ZPcisSNP(y, X=cbind(1, X, g), G=G, PXEM=TRUE, maxIter=1000)
simLike <- LRTSimZP(Z = X, E = g, G = G, nsim=1e5, parallel=c("multicore"), ncpus = 4L) ## exact LRT

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
Haojie Lu, and Ping Zeng<sup>#</sup> (2021). Leveraging trans-ethnic genetic risk scores to improve association power for complex traits in underrepresented populations.

# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn.

# Update
2022-07-22 GAMM version 1.0
