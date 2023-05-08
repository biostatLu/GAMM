// #include <RcppArmadillo.h>
#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>

using namespace Rcpp;
using namespace arma;
using namespace std;

//' @param y  gene expression vector
//' @param X  covariates file for eQTL data
//' @param G  normalized genotype (cis-SNPs) matrix for eQTL
//' @param maxIter  maximum iteration (default is 1000)

void lmm_pxem_ptr2_ZPcisSNP(const arma::vec& y, 
                      const arma::mat& X, 
                      const arma::mat& G, 
                      const int& maxIter,
                      double& sigma2e,
                      double& sigma2b,
                      arma::vec& betax,
                      double& loglik_max,
                      int& iteration,
					  arma::mat& covb,
					  arma::vec& mub,
					  bool PXEM) {

  int n = y.n_elem, p = G.n_cols;

  if (y.n_elem != G.n_rows || G.n_rows != X.n_rows){
     perror("The dimensions in outcome and covariates (X and G) are not matched");
 }

  /*
  if (betax.n_elem != X.n_cols){
    perror("The dimensions in covariates are not matched in X and betax");
  }

  if (p != mub.n_elem){
    perror("The dimensions in covariates are not matched in mub");
  }

  if (p != Sigb.n_cols){
    perror("The dimensions in covariates are not matched in Sigb");
  }
  */

  mat GtG = G.t()*G, XtX = X.t()*X, XtG = X.t()*G;
  vec Gty = G.t()*y, Xty = X.t()*y;

  vec SXy;
  mat SXG;

  if(X.n_cols==1){
    SXy = mean(y);
    SXG = mean(G,0);
  } 
  else{
    SXy = solve(XtX, Xty);
    SXG = solve(XtX, XtG);
  }

  double gam, gam2;  // parameter expansion
  vec eVal;
  mat eVec;

  eig_sym(eVal, eVec, GtG);

  // initialize
  sigma2e = var(y);
  sigma2b = sigma2e/p;
  betax = SXy - SXG * mub;
  vec loglik(maxIter);
  loglik(0) = -datum::inf;

  vec D;
  vec Gmu;
  vec y_bar = y - X * betax;
  double y_Gmu2, E, pix = n/2 * log(2 * datum::pi);

  iteration = maxIter-1;

  //PXEM algorithm
  if (PXEM) {
  for (int iter = 1; iter < maxIter; iter ++ ) {
    // E-step
      D = 1/sigma2b + eVal/sigma2e;
    mub = 1/sigma2e * eVec * (eVec.t() * (G.t() * y_bar)/D);
    Gmu = G * mub;
    y_Gmu2 = sum(pow(y_bar - Gmu, 2));

    // Evaluate loglik
    E = y_Gmu2/(2 * sigma2e) + accu(pow(mub, 2))/(2 * sigma2b);
    loglik(iter) = - p * log(sigma2b)/2 - n * log(sigma2e)/2 - E - sum(log(D))/2 - pix;

    if ( loglik(iter) - loglik(iter - 1) < 0 ){
      perror("The likelihood failed to increase!");
    }

    if (abs(loglik(iter) - loglik(iter - 1)) < 1e-10) {
      iteration = iter;
      break;
    }

    // M-step
    gam = sum(y_bar % Gmu) / (accu(pow(Gmu, 2)) + sum(eVal/D));
    gam2 = pow(gam, 2);

    betax = SXy - (SXG * mub) * gam;
    y_bar = y - X * betax;;

    sigma2e = sum(pow(y_bar - Gmu * gam, 2))/n + gam2 * sum(eVal/D)/n;
    sigma2b = accu(pow(mub, 2))/p + sum(1/D)/p;

    // Reduction step
    sigma2b = gam2 * sigma2b;
    // gam = 1;
    // gam2 = pow(gam , 2);
  }}


    //standard EM algorithm
  else {
   for (int iter = 1; iter < maxIter; iter ++ ) {
    // E-step
      D = 1/sigma2b + eVal/sigma2e;
    mub = 1/sigma2e * eVec * (eVec.t() * (G.t() * y_bar)/D);
    Gmu = G * mub;
    y_Gmu2 = sum(pow(y_bar - Gmu, 2));

    // Evaluate loglik
    E = y_Gmu2/(2 * sigma2e) + accu(pow(mub, 2))/(2 * sigma2b);
    loglik(iter) = - p * log(sigma2b)/2 - n * log(sigma2e)/2 - E - sum(log(D))/2 - pix;

    if ( loglik(iter) - loglik(iter - 1) < 0 ){
      perror("The likelihood failed to increase!");
    }

    if (abs(loglik(iter) - loglik(iter - 1)) < 1e-10) {
      iteration = iter;
      break;
    }

    // M-step
    //gam = sum(y_bar % Gmu) / (accu(pow(Gmu, 2)) + sum(eVal/D));
    //gam2 = pow(gam, 2);
    betax = SXy - (SXG * mub);
    y_bar = y - X * betax;
    sigma2e = sum(pow(y_bar - Gmu, 2))/n + sum(eVal/D)/n;
    sigma2b = accu(pow(mub, 2))/p + sum(1/D)/p;
    // Reduction step
    //sigma2b = gam2 * sigma2b;
    // gam = 1;
    // gam2 = pow(gam , 2);
   }}

  vec loglik_out;
  loglik_out = loglik.subvec(0, iteration);
  loglik_max = loglik(iteration);

  ///////////////////
  mat Ip = eye(p, p);
  if (sigma2b < 1E-10) {covb = sigma2e * inv(XtX);}
  else {covb = sigma2e * inv(XtX - XtG * inv(Ip * sigma2e/sigma2b + GtG) * XtG.t());}

			//for(size_t i=0; i<2; ++i) {
			//	for(size_t j=0; j<2; ++j) {
			//		cout<<covb(i,j)<<"\t";
			//	}
			//	cout<<endl;	
			//}
}

//' @param y  gene expression vector.
//' @param X  covariates file for eQTL data.                  
//' @param G  normalized genotype (cis-SNPs) matrix for eQTL.    
//' @param maxIter  maximum iteration (default is 1000).


// [[Rcpp::export]]
SEXP lmm_pxem2_ZPcisSNP(const arma::vec y, const arma::mat X, const arma::mat G, const bool PXEM, const int maxIter){

    double sigma2e = var(y)/2, sigma2b = var(y)/2, loglik;
    vec betax =zeros<vec>(X.n_cols); 
    int iter;
	mat covb = zeros<mat>(X.n_cols, X.n_cols);
    vec mub  = zeros<vec>(G.n_cols);

    lmm_pxem_ptr2_ZPcisSNP(y, X, G, maxIter, sigma2e, sigma2b, betax, loglik, iter, covb, mub, PXEM);
	vec theta = zeros<vec>(2);
	theta(0) = sigma2e;
	theta(1) = sigma2b;

    List output = List::create(Rcpp::Named("sigma2e")    = sigma2e,
                               Rcpp::Named("sigma2b")    = sigma2b,
                               Rcpp::Named("theta")      = theta,
		                       Rcpp::Named("alpha")      = betax,
                               Rcpp::Named("loglik")     = loglik,
                               Rcpp::Named("iteration")  = iter,
                               Rcpp::Named("cov")        = covb,
							   Rcpp::Named("mub")        = mub);
    return output;
}
