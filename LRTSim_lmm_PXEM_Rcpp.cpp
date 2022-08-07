
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]
#include <fstream>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>

using namespace std;
using namespace arma;
using namespace Rcpp;

#define ARMA_DONT_PRINT_ERRORS

// [[Rcpp::export]]
SEXP PX(SEXP P0, SEXP sqrtK)//这里的特征值按照从小到大排序，与eigen相反
{
	try {
		mat K = as<mat>(sqrtK);
		mat P = as<mat>(P0);
		K = K * P * K;
		return List::create(Named("PX") = K);

	} catch (std::exception &ex) {
		forward_exception_to_r( ex );
	} catch (...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}

// [[Rcpp::export]]
SEXP eigenK(SEXP Kmat, bool sqrtk)//这里的特征值按照从小到大排序，与eigen相反
{
	try {
		mat K = as<mat>(Kmat);
		mat vectors;
		vec eigval;

		eig_sym(eigval, vectors, K, "dc");
		const size_t num = sum(eigval < 1e-5);
		const uvec  idx = find(eigval < 1e-5);
		for (size_t i=0; i<num; i++) {eigval[idx[i]] = 0;}

		if (sqrtk) {
		mat sqrtK = vectors * diagmat(sqrt(eigval)) * vectors.t();
		return List::create(Named("value") = eigval, Named("vectors") = vectors, Named("sqrtK") = sqrtK);
		}

		else {
		return List::create(Named("value") = eigval, Named("vectors") = vectors);
		}

	} catch (std::exception &ex) {
		forward_exception_to_r( ex );
	} catch (...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}


// [[Rcpp::export]]
List RLRsimCpp (
	int p,
	int k,
	int n,
	int nsim,
	int g,
	int q,
	Rcpp::NumericVector mu,
	Rcpp::NumericVector lambda,
	double lambda0,
	Rcpp::NumericVector xi,
	bool REML) {

	Rcpp::RNGScope scope;

	/* allocate: */
	Rcpp::NumericMatrix lambdamu(g, k);
	Rcpp::NumericMatrix lambdamuP1(g, k);
	Rcpp::NumericMatrix fN(g, k);
	Rcpp::NumericMatrix fD(g, k);

	Rcpp::NumericVector sumlog1plambdaxi(g);
	Rcpp::NumericVector Chi1(k);
	Rcpp::NumericVector res(nsim);
	Rcpp::IntegerVector lambdaind(nsim);

	int is, ig, ik, dfChiK, n0, kx;
	double  LR, N, D, ChiK, ChiSum;

	kx = k;
	if (k == n) {kx = n - p;}
	//cout << k << " " << n << " " << kx << endl;
	dfChiK = n-p-k; if (dfChiK < 0) { dfChiK = 0;};
	if (REML) {
	n0 = n - p;
	for (ik=0; ik < k; ++ik) {
		xi[ik] = mu[ik];
		}
	}	else {
	n0 = n;
	}

	/*precompute stuff that stays constant over simulations*/
	for (ig = 0; ig < g; ++ig) {
		sumlog1plambdaxi[ig] = 0;
		for (ik=0 ; ik < k; ++ik) {
			lambdamu(ig, ik) = lambda[ig] * mu[ik];
			lambdamuP1(ig, ik) = lambdamu(ig, ik) + 1.0;

			//fN(ig, ik) = ((lambda[ig] - lambda0) * mu[ik]) / lambdamuP1(ig, ik);
			//fD(ig, ik) = (1 + lambda0 * mu[ik]) / lambdamuP1(ig, ik);

			fN(ig, ik) = (lambda[ig] * mu[ik]) / lambdamuP1(ig, ik);
			fD(ig, ik) = 1 / lambdamuP1(ig, ik);

			sumlog1plambdaxi[ig] += log1p(lambda[ig] * xi[ik]);
		}/* end for k*/
	}/* end for g*/

	for (is = 0; is < nsim; ++is) {
	/*make random variates, set LR 0*/
		LR =  0;
		//ChiSum = 0;
		ChiK = rchisq(1, dfChiK)[0];//[0]必须有,等价于rchisq(1, dfChiK)
		//cout << ChiK << " " << dfChiK << endl;
		//dfChiK = n-p-k; if (dfChiK < 0){ dfChiK = 0;};
		Chi1 = rchisq(k, 1);
		//if (!REML) {ChiSum = std::accumulate(Chi1.begin(), Chi1.end(), 0.0);}//等价于sum(Chi1)
		//ChiSum = std::accumulate(Chi1.begin(), Chi1.end(), 0.0);
		for (ig = 0; ig < g; ++ig) {
		/*loop over lambda-grid*/
			N = D = 0;
			//for (ik=0; ik < k; ++ik) {
			for (ik=0; ik < kx; ++ik) {
				/*loop over mu, xi*/
				N = N + fN(ig, ik) * Chi1[ik];
				D = D + fD(ig, ik) * Chi1[ik];
				}
			D = D + ChiK;
			LR = n0 * log1p(N/D) - sumlog1plambdaxi[ig];
			if (LR >= res[is]) {
			/*save if LR is bigger than previous LR*/
				res[is] = LR;
				lambdaind[is] = ig + 1;
				} else break;
		}/*end for g*/
		/* add additional term for LR*/
		//cout << res[is] << endl;
		//if (!REML){res[is] = res[is] + n * log1p(rchisq(1, q)[0] / (ChiSum + ChiK));}
		//if (!REML){res[is] = res[is] + n * log1p(rchisq(1, q)[0] / rchisq(1, n-p)[0]);}
		res[is] = res[is] + n * log1p(rchisq(1, q)[0] / rchisq(1, n-p)[0]);
		//res[is] = res[is] + n * log1p(rchisq(1, q)[0] / (ChiSum + ChiK));
		//cout << res[is] << endl;
		//if (is % (nsim/5) == 0) {cout << (is/nsim)*100 << " -- finished-- " << endl;}
		}/*end for nsim*/

	return List::create(Named("res") = res,
	Named("lambdaind") = lambdaind,
	Named("lambdamu") = lambdamu,
	Named("fN") = fN,
	Named("fD") = fD,
	Named("sumlog1plambdaxi") = sumlog1plambdaxi,
	Named("Chi1") = Chi1,
	Named("ChiK") = ChiK);
}
