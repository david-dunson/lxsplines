#include <RcppArmadillo.h>
#include <iostream>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
extern "C" {
	//////////////////////////////////////////////////////////////////
	// Stub for 
	// Bivariate and Trivariate Normal integration routines 
	// currently in the Fortran Files
	double  bvnd_(double * DH, double *DK, double *R);
	double  tvtl_(int * NU, double *HK, double *R, double *EPSI);
	double  sadmvn_(int * N, double * LOWER, double* UPPER,
		int * INFIN, double* CORREL, int *MAXPTS,
		double * ABSEPS, double * RELEPS, double * E, double * VALUE, int *INFORM);
	//////////////////////////////////////////////////////////////////
}


	/**********************************************************
	*Funtion BivNProb
	*Input: double *mean -- mean vector (implicilty vector of 3)
	*		double *cv   -- covariance in tri-diagonal form First i.e.
	1  x
	2  3
	where numbers represent the array index
	double *RV   -- pointer to a return value
	*Return: The probability that X > 0
	*Note:    Assumes all the numbers passed work, no error checking.
	**************************************************************/
	int BivNProb(double *mean, double *cv, double *RV) {
		double a = -1.0*mean[0] / sqrt(cv[0]);
		double b = -1.0*mean[1] / sqrt(cv[3]);
		double rho = cv[1] / sqrt(cv[0] * cv[3]);
		*RV = bvnd_(&a, &b, &rho);
		return 0;
	}

	

	/******************************************************************
	*Funtion TriNProb
	*Input: double *mean -- mean vector (implicilty vector of 3)
	*		double *cv   -- covariance in tri-diagonal form First i.e
	1  x  x
	2  4  x
	3  5  6
	where numbers represent the array index
	double *RV   -- pointer to a return value
	*Return: The probability that X > 0
	*Note:    Assumes all the numbers passed work, no error checking.
	************************************************************/
	int TriNProb(double *mean, double *cv, double *RV) {
		double m[3];
		double R[3];
		int df = 0;
		double EPS = 1e-14;

		//  compute the means
		m[0] = double(mean[0]) / sqrt(cv[0]);
		m[1] = double(mean[1]) / sqrt(cv[4]);
		m[2] = double(mean[2]) / sqrt(cv[8]);

		// compute the correlations
		R[0] = cv[1] / sqrt(cv[0] * cv[4]);
		R[1] = cv[2] / sqrt(cv[0] * cv[8]);
		R[2] = cv[5] / sqrt(cv[4] * cv[8]);

		*RV = tvtl_(&df, m, R, &EPS);
		return 0;
	}

	

	/*********************************************************
	Function:genTruncNormZ
	Purpose:  Generate a truncated normal rv, where the truncation is at zero
	Input  :  mean - Mean value
			  sd   - standard deviation
			  rV   - pointer to the return value
	output : rV    - a truncated normal random variable 	
	*********************************************************/
	int genTruncNormZ(double* mean, double* sd, double* rV) {
		if (*mean < 0) {
			// do everything on the log scale  - more numerically stable
			double ucutoff = -1.0*R::pnorm(0.0, -1.0*(*mean), *sd, true, true);
			GetRNGstate();
			// crazy numerics go like this: 
			// -log(uniform) = exponential distribution
			// exponential is memoryless so the cuttoff can be moved using the following
			// easy peasy code
			double cv = R::rexp(1) + ucutoff;
			// once we have the new cut off we back translate
			PutRNGstate();
			*rV = -1.0*R::qnorm(-1.0*cv, -1.0*(*mean), *sd, true, true);
			if (*rV < 0) { *rV = 1e-5; }//very rare numerical problem
		}
		else {
			// no real problems when 0 < mean
			// do it the good old fasion way
			double lcutoff = R::pnorm(0.0, *mean, *sd, true, false);
			GetRNGstate();
			double tunif = R::runif(lcutoff, 1.0);
			PutRNGstate();
			*rV = R::qnorm(tunif, *mean, *sd, true, false);
		}

		return  0;
	}



// function that calls Fortran multivariate normal library
// for 4-variate multivariate normal
// [[Rcpp::export]]
double SADMVN(arma::mat M, arma::mat C) {
	int    N = 4;
	int    MAXPTS = 2000;
	double RELEPS = 0;
	double ABSEPS = 1e-6;
	double E_E = 0.0;
	double returnV = 0.0;
	double UPPER[4] = { std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity() }; //POPULATED NEVER USED
	int INFIN[4] = { 1,1,1,1 };
	int INFORM;


	// Set up the correlation matrix
	int gk = 0;
	double B[6];
	for (int i = 1; i < 4; i++) {
		M(i - 1, 0) = M(i - 1, 0) / sqrt(C(i - 1, i - 1));
		for (int j = 0; j < i; j++) {
			B[gk] = C(i, j) / sqrt(C(j, j)*C(i, i));
			gk++;
		}
	}
	M(3, 0) = M(3, 0) / sqrt(C(3, 3));

	M = -M;
	sadmvn_(&N, M.memptr(), UPPER, INFIN, B, &MAXPTS, &ABSEPS, &RELEPS, &E_E, &returnV, &INFORM);

	return returnV;


}

double rtn(double mean, double sd) {
	double rV = 0.0;

	if (mean < 0) {
		// do everything on the log scale  - more numerically stable
		double ucutoff = -1.0*R::pnorm(0.0, -1.0*(mean), sd, true, true);
		GetRNGstate();
		double cv = R::rexp(1) + ucutoff;
		PutRNGstate();
		rV = -1.0*R::qnorm(-1.0*cv, -1.0*(mean), sd, true, true);
	}
	else {
		double lcutoff = R::pnorm(0.0, mean, sd, true, false);
		GetRNGstate();
		double tunif = R::runif(lcutoff, 1.0);
		PutRNGstate();
		rV = R::qnorm(tunif, mean, sd, true, false);
	}

	return  rV;
}

// [[Rcpp::export]]
double TriNProb(NumericVector mean, NumericVector cv) {
	// Function BivNProb:
	//Wrapper function that computes the probability  that X > 0
	// by calling the fortran function bvnd_
	// assumes everythign is great no error checking
	double m[3];
	double R[3];
	double RV;
	int df = 0;
	double EPS = 2e-16;

	//  compute the means
	m[0] = double(mean[0]) / sqrt(cv[0]);
	m[1] = double(mean[1]) / sqrt(cv[4]);
	m[2] = double(mean[2]) / sqrt(cv[8]);

	// compute the correlations
	R[0] = cv[1] / sqrt(cv[0] * cv[4]);
	R[1] = cv[2] / sqrt(cv[0] * cv[8]);
	R[2] = cv[5] / sqrt(cv[4] * cv[8]);

	RV = tvtl_(&df, m, R, &EPS);
	return RV;
}

// [[Rcpp::export]]
double BivNProb(NumericVector mean, NumericVector cv) {
	// Function BivNProb:
	//Wrapper function that computes the probability  that X > 0
	// by calling the fortran function bvnd_
	// assumes everythign is great no error checking
	double rv;
	double a = -1.0*mean[0] / sqrt(cv[0]);
	double b = -1.0*mean[1] / sqrt(cv[3]);
	double rho = cv[1] / sqrt(cv[0] * cv[3]);
	rv = bvnd_(&a, &b, &rho);
	return rv;
}

// [[Rcpp::export]]
NumericVector rtmvn(NumericVector tMean, NumericVector tVar) {
	arma::mat COV = Rcpp::as<arma::mat>(tVar);
	arma::mat MEAN = Rcpp::as<arma::mat>(tMean);
	arma::mat tV12_22inv_tM2(1, MEAN.n_rows);
	arma::mat tV(1, MEAN.n_rows);
	arma::mat tV12_22inv(1, MEAN.n_rows*(MEAN.n_rows - 1));
	arma::mat RV = MEAN;



	if (MEAN.n_rows == 1) {
		RV(0, 0) = rtn(MEAN(0, 0), sqrt(COV(0, 0)));
	}
	else {
		for (int i = 0; i < int(MEAN.n_rows); i++) {
			arma::mat TEMP = COV;
			arma::mat TEMP2 = COV.cols(i, i);
			arma::mat TEMP3 = MEAN;

			TEMP3.shed_row(i);
			TEMP.shed_row(i);
			TEMP.shed_col(i);

			TEMP2.shed_row(i);
			tV12_22inv.submat(0, i*(MEAN.n_rows - 1), 0, (i + 1)*(MEAN.n_rows - 1) - 1) = TEMP2.t()*TEMP.i();
			tV.submat(0, i, 0, i) = sqrt(COV.submat(i, i, i, i) - tV12_22inv.submat(0, i*(MEAN.n_rows - 1), 0, (i + 1)*(MEAN.n_rows - 1) - 1)*TEMP2); // variances conditional on the mean
			tV12_22inv_tM2.submat(0, i, 0, i) = MEAN.submat(i, 0, i, 0) - TEMP2.t()*solve(TEMP, TEMP3);
		}
		// set the starting value

		for (int i = 0; i< int(MEAN.n_rows); i++) {
			if (RV(i, 0) < 0.0) {
				RV(i, 0) = 0.01;
			}
		}

		arma::mat TempV = RV;


		for (int i = 0; i < 20; i++) {
			// gibbs sampler
			for (int j = 0; j < int(MEAN.n_rows); j++) {
				TempV = RV;
				TempV.shed_row(j);
				TempV = tV12_22inv_tM2.submat(0, j, 0, j) + tV12_22inv.submat(0, j*(MEAN.n_rows - 1), 0, (j + 1)*(MEAN.n_rows - 1) - 1)*TempV;
				RV(j, 0) = rtn(TempV(0, 0), tV(0, j));
			}
		}

	}


	// return the R value
	return wrap(RV);
}

// [[Rcpp::export]]
SEXP qcopy(SEXP a, SEXP b, SEXP c, SEXP d) {
	double *x = REAL(a);
	double *y = REAL(b);
	int n = INTEGER(c)[0];
	int i = INTEGER(d)[0];
	memmove(&x[n*(i - 1)], y, n * sizeof(double));
	return R_NilValue;
}
