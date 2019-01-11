//////////////////////////////////////////////////////////////////////
//
//
//
//
//
//
//////////////////////////////////////////////////////////////////////
#include <RcppArmadillo.h>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;


double SADMVN(arma::mat M, arma::mat C); 
double TriNProb(NumericVector mean, NumericVector cv); 
double BivNProb(NumericVector mean, NumericVector cv); 

int BivNProb(double *mean, double *cv, double *RV);
int TriNProb(double *mean, double *cv, double *RV);

int genTruncNormZ(double* mean, double* sd, double* rV);
double rtn(double mean, double sd); 

NumericVector rtmvn(NumericVector tMean, NumericVector tVar);

/******************************************************************
function: sampleBetas()
purpose:  Gibb Sampler of the betas one at a time
Input  :  ttY - Data Vector
		  ttX - X covariate vector
		  tbetas - current betas
		  omega  - for the weighted regression N x 1
Return :  Gibbs Sample of the Betas 
*******************************************************************/
// [[Rcpp::export]]
NumericVector sampleBetas(NumericVector ttY, NumericVector ttX, NumericVector tbetas,
						  NumericVector LAM, NumericVector intLAM, NumericVector p,
						  NumericVector tau, NumericVector omega) {

	arma::mat Y = Rcpp::as<arma::mat>(ttY);
	arma::mat X = Rcpp::as<arma::mat>(ttX);
	arma::mat W = Rcpp::as<arma::mat>(omega); 
	arma::mat betas = Rcpp::as<arma::mat>(tbetas);
	
	//////////////////////////////////////////////////////////////
	// first beta is the background
	arma::mat tX = X.col(0); // first column
	arma::mat tempX = tX; 
	tempX.each_col() %=  W;

	Y = Y - X.submat(0, 1, X.n_rows - 1, X.n_cols - 1) * betas.submat(1, 0, X.n_cols - 1, 0);
	
	arma::mat tV = 1 / ((tau[0])*(tempX.t()*tX) + intLAM[0]);

	arma::mat tM = ((tau[0])*(tempX.t()*Y) + intLAM[0] * 0.0)*tV;
	
	GetRNGstate();
	betas(0, 0) = R::rnorm(tM(0, 0), sqrt(tV(0, 0)));
	PutRNGstate();
	
	Y = Y - tX * betas(0, 0);
	double PL = LAM[0];
	/////////////////////////////////////////////////////////
	//  Sample  each  constrained beta one at a time
	for (signed int i = 1; i < X.n_cols; i++) {

		tX = X.submat(0, i, X.n_rows - 1, i);
		tempX = tX ; 	tempX.each_col() %= W;
		Y = Y + tX * betas(i, 0); // the  previous iteration this value was r
								  // removed from  the  y vector  'add it back'
		tV = ((tau[0])*(tempX.t()*tX)).i();
		tM = ((tau[0])*(tempX.t()*Y) - PL)*tV;

		double A = log(p[0]) + R::dnorm(0.0, tM(0, 0), sqrt(tV(0, 0)), true);
		double B = log((1.0 - p[0])*(PL)) + R::pnorm(0.0, tM(0, 0), sqrt(tV(0, 0)), false, true);

		if (A > B) {
			A = A - A;
			B = B - A;
		}else {
			A = A - B;
			B = B - B;
		}

		double PZERO = exp(A) / (exp(A) + exp(B));
		
		if (R::runif(0.0, 1.0) < PZERO) {
			betas(i, 0) = 0.0;
		}
		else {
			double rV = 0.0;
			double mean = tM(0, 0); double sd = sqrt(tV(0, 0));
			genTruncNormZ(&mean, &sd, &rV);
			if (rV < 0.0) { rV = 1e-16; }
			betas(i, 0) = rV;
		}
		Y = Y - tX * betas(i, 0); // add the new *regressor* to the residual
	}
	return wrap(betas);
}

/***********************************************************************
function: sinsertBeta_trivariate 
Purpose : Gibbs Sampler for a Knot Insert of the Betas- Trivariate Algorithm
Return  : The values of the Betas 
************************************************************************/
// [[Rcpp::export]]
List sinsertBeta_trivariate(NumericVector tY, NumericVector Xidx, NumericVector ctau,
							NumericVector tp, NumericVector lam,  NumericVector omega) {

	// initialize the proper values
	double p = tp[0];
	arma::mat Y = Rcpp::as<arma::mat>(tY);
	arma::mat X = Rcpp::as<arma::mat>(Xidx);
	arma::mat W = Rcpp::as<arma::mat>(omega);
	double tau = ctau[0];

	arma::mat tLAM = Rcpp::as<arma::mat>(lam);
	arma::mat LAM(tLAM.n_rows, tLAM.n_rows);  LAM.diag() = tLAM;
	arma::mat INTC(16, 2);
	arma::mat PMAT(16, 1);

	INTC.zeros();
	arma::mat XW = X; 
	XW.each_col() %= W;

	arma::mat B = (tau)*(XW.t()*X);          arma::mat U = (tau)*(XW.t()*Y) - lam[0];

	List uniMeans(4); List uniVars(4);
	List bivMeans(6); List bivVars(6);
	List TriMeans(4); List TriVars(4);

	List returnV(4);
	StringVector names(4);

	names(0) = "VAR";    names(1) = "MEAN";
	names(2) = "LPROB";  names(3) = "S_ELM";

	returnV.attr("names") = names;

	double cnum = arma::cond(B);

	if (cnum > 1e5) {
		returnV[0] = wrap(0);
		returnV[1] = wrap(0);
		returnV[2] = wrap(log(0));
		returnV[3] = wrap(1);
		return returnV;
	}

	PMAT(0, 0) = (p)*(p)*(p)*(p);

	arma::mat tV;            arma::mat tM;
	arma::mat TEMP; double pr = 0.0;
	// The  univariate values
	for (int i = 0; i < 4; i++) {
		PMAT(i + 1, 0) = p * p*p*(1 - p);
		tV = (B(i, i));
		tV = tV.i();
		tM = tV * (U(i, 0));
		TEMP = -.5*tM.t()*solve(tV, tM);
		INTC(i + 1, 0) = TEMP(0, 0);
		//Total positive area for the univariate normal distribution
		pr = 1.0 - R::pnorm(0, tM(0, 0), sqrt(tV(0, 0)), true, false);

		if (pr <= 0 || isnan(pr)) { pr = 0; }// rare numerical error when probability is essentially zero
		TEMP = log(pr*sqrt(tV)*lam[0]) + 0.5*log(2.0*M_PI);
		INTC(i + 1, 1) = TEMP(0, 0);
		uniMeans[i] = wrap(tM(0, 0)); uniVars[i] = wrap(tV(0, 0));
	}
	//bivariate normal case
	arma::mat tU;
	arma::mat tB;
	int counter = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 4; j++) {
			PMAT(counter + 5, 0) = p * p*(1 - p)*(1 - p);
			tU = join_cols(U.row(i), U.row(j));
			tB = join_rows(B.col(i), B.col(j));
			for (int k = 3; k >= 0; k--) {
				if (k != i && k != j)
					tB.shed_row(k);
			}
			tV = (tB);     tV = tV.i();
			tM = tV * (tU);

			bivVars[counter] = wrap(tV);
			bivMeans[counter] = wrap(tM);
			TEMP = -.5*tM.t()*solve(tV, tM);
			INTC(counter + 5, 0) = TEMP(0, 0);
			pr = BivNProb(wrap(tM), wrap(tV));  if (pr <= 0 || isnan(pr)) { pr = 0; }// numerical instability for small probabilities
			TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0] * lam[0])) + log(2.0*M_PI);

			INTC(counter + 5, 1) = TEMP(0, 0);
			counter++;
		}
	}

	//trivariate case 
	for (int i = 0; i<4; i++) {
		PMAT(i + 11, 0) = p * (1 - p)*(1 - p)*(1 - p);
		tU = U; tU.shed_row(i);
		tB = B; tB.shed_col(i); tB.shed_row(i);
		tV = (tB);    tV = tV.i();
		tM = tV * (tU);
		TriVars[i] = wrap(tV);
		TriMeans[i] = wrap(tM);
		TEMP = -.5*tM.t()*solve(tV, tM);
		INTC(i + 11, 0) = TEMP(0, 0);
		pr = TriNProb(wrap(tM), wrap(tV));   if (pr <= 0 || isnan(pr)) { pr = 0; }// numerical instability for small probabilities

		TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0] * lam[0] * lam[0])) + 1.5*log(2.0*M_PI);
		INTC(i + 11, 1) = TEMP(0, 0);
	}

	// 4-variate normal
	PMAT(15, 0) = (1 - p)*(1 - p)*(1 - p)*(1 - p);
	tU = U;    tB = B;
	tV = (tB);    tV = tV.i();
	tM = tV * (tU);
	TEMP = -.5*tM.t()*solve(tV, tM);
	INTC(15, 0) = TEMP(0, 0);
	pr = SADMVN(tM, tV);   if (pr <= 2e-16 || isnan(pr)) { pr = 0; } // rare numerical error caused by low probabilities and a large 
																	 // unstable determinant 

	TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0] * lam[0] * lam[0] * lam[0])) + 2 * log(2.0*M_PI);
	if (TEMP(0, 0) > 10) { TEMP(0, 0) = log(0); } // another rare numerical problem
	INTC(15, 1) = TEMP(0, 0);

	arma::mat TT = INTC.col(1) - INTC.col(0);

 
	arma::mat MT = max(TT);
	arma::mat SMULT = PMAT % exp(TT - repmat(MT, 16, 1));
	arma::mat log_prob = -0.5*(tau)*(Y.t()*Y) + log(sum(SMULT)) + MT;
	arma::mat POST_PROBS = exp(log(SMULT) - log(repmat(sum(SMULT), 16, 1))); //posterior sampling prob
	GetRNGstate();
	// get the cumulative  density
	double rvalue = R::runif(0, 1); int rval = -1;

	for (int i = 1; i < 16; i++) { POST_PROBS(i, 0) += POST_PROBS(i - 1, 0); }
	for (int i = 0; i < 16 && rval == -1; i++) { if (rvalue <= POST_PROBS(i, 0)) { rval = i; } }


	rval = rval + 1;
	PutRNGstate();
	NumericVector rVar, rMeans;

	switch (rval) {
	case 16:
		rVar = wrap(tV); rMeans = wrap(tM);
		break;
	case 1:
		rVar = 0.0;    rMeans = 0.0;
		break;
	case 2: case 3: case 4: case 5: // univariate chosen
		rVar = uniVars[rval - 2];
		rMeans = uniMeans[rval - 2];
		break;
	case 12: case 13: case 14: case 15:
		rVar = TriVars[rval - 12];
		rMeans = TriMeans[rval - 12];
		break;
	default:
		rVar = bivVars[rval - 6];
		rMeans = bivMeans[rval - 6];
		break;
	}


	returnV[0] = wrap(rVar);
	returnV[1] = wrap(rMeans);
	returnV[2] = wrap(log_prob(0, 0));
	returnV[3] = wrap(rval);



	return returnV;
}


/***********************************************************************
function: sdeleteBeta_trivariate
Purpose : Gibbs Sampler for a Knot Delete of the Betas- Trivariate Algorithm
Return  : The values of the Betas
************************************************************************/
// [[Rcpp::export]]
List sdeleteBeta_trivariate(NumericVector tY, NumericVector Xidx, NumericVector ctau,
	                          NumericVector tp, NumericVector lam,  NumericVector omega) {

	// initialize the proper values
	arma::mat Y = Rcpp::as<arma::mat>(tY);
	arma::mat X = Rcpp::as<arma::mat>(Xidx);
	arma::mat W = Rcpp::as<arma::mat>(omega);
	double tau = ctau[0];

	arma::mat tLAM = Rcpp::as<arma::mat>(lam);
	arma::mat LAM(tLAM.n_rows, tLAM.n_rows);  LAM.diag() = tLAM;
	arma::mat INTC(8, 2);

	INTC.zeros();
	arma::mat XW = X;
	XW.each_col() %= W;

	arma::mat B = (tau)*(XW.t()*X);          arma::mat U = (tau)*(XW.t()*Y) - lam[0];
	double p = tp[0];

	arma::mat tV;            arma::mat tM;
	arma::mat uniMeans(3, 1); arma::mat uniVars(3, 1);
	arma::mat bivMeans(2, 3); arma::mat bivVars(2, 6);

	arma::mat TEMP; double pr = 0.0;
	arma::mat PMAT(8, 1); PMAT(0, 0) = (p)*(p)*(p);

	// start with the univariate values
	for (int i = 1; i <4; i++) {
		PMAT(i, 0) = p * p*(1 - p);
		tV = (B(i - 1, i - 1));
		tV = tV.i();
		tM = tV * (U(i - 1, 0) - lam[0]);
		TEMP = -.5*tM.t()*solve(tV, tM);
		INTC(i, 0) = TEMP(0, 0);
		//Total positive area for the univariate normal distribution
		pr = 1 - R::pnorm(0, tM(0, 0), sqrt(tV(0, 0)), true, false);
		TEMP = log(pr*sqrt(tV)*lam[0]) + 0.5*log(2.0*M_PI);
		INTC(i, 1) = TEMP(0, 0);
		uniMeans(i - 1, 0) = tM(0, 0); uniVars(i - 1, 0) = tV(0, 0);

	}
	//  bivariate values don't worry about looping as it just adds extra
	//  First  Condition


	PMAT(4, 0) = p * (1 - p)*(1 - p);
	arma::mat tU = U.submat(0, 0, 1, 0);
	arma::mat tB = B.submat(0, 0, 1, 1);
	tLAM = LAM.submat(0, 0, 1, 1);
	tV = (tB);
	tV = tV.i();
	tM = tV * (tU - lam[0]);
	bivVars.submat(0, 0, 1, 1) = tV;
	bivMeans.submat(0, 0, 1, 0) = tM;
	TEMP = -.5*tM.t()*solve(tV, tM);
	INTC(4, 0) = TEMP(0, 0);
	pr = BivNProb(wrap(tM), wrap(tV));
	if (pr < 0 || isnan(pr)) { pr = 0.0; } // rare numerical error
	TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0] * lam[0])) + log(2.0*M_PI);
	INTC(4, 1) = TEMP(0, 0);

	// Second Condition

	PMAT(5, 0) = p * (1 - p)*(1 - p);
	tB(0, 0) = B(0, 0); tB(0, 1) = B(0, 2);
	tB(1, 0) = B(2, 0); tB(1, 1) = B(2, 2);
	tLAM(0, 0) = lam[0]; tLAM(1, 1) = lam[2];
	tU(0, 0) = U(0, 0); tU(1, 0) = U(2, 0);
	tV = (tB);
	tV = tV.i();
	tM = tV * (tU - lam[0]);
	bivVars.submat(0, 2, 1, 3) = tV;
	bivMeans.submat(0, 1, 1, 1) = tM;
	TEMP = -.5*tM.t()*solve(tV, tM);
	INTC(5, 0) = TEMP(0, 0);

	pr = BivNProb(wrap(tM), wrap(tV));
	if (pr < 0 || isnan(pr)) { pr = 0.0; } // rare numerical error
	TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0] * lam[0])) + log(2.0*M_PI);
	INTC(5, 1) = TEMP(0, 0);

	// Third Condition
	PMAT(6, 0) = p * (1 - p)*(1 - p);
	tU = U.submat(1, 0, 2, 0);
	tB = B.submat(1, 1, 2, 2); tLAM = LAM.submat(1, 1, 2, 2);
	tV = (tB);
	tV = tV.i();
	tM = tV * (tU - lam[0]);
	bivVars.submat(0, 4, 1, 5) = tV;
	bivMeans.submat(0, 2, 1, 2) = tM;
	TEMP = -.5*tM.t()*solve(tV, tM);
	INTC(6, 0) = TEMP(0, 0);
	pr = BivNProb(wrap(tM), wrap(tV));
	if (pr < 0 || isnan(pr)) { pr = 0.0; } // rare numerical error
	TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0] * lam[0])) + log(2.0*M_PI);
	INTC(6, 1) = TEMP(0, 0);

	// TRIVARIATE CASE
	PMAT(7, 0) = (1 - p)*(1 - p)*(1 - p);
	tV = (B);
	tV = tV.i();
	tM = tV * (U - lam[0]);
	TEMP = -.5*tM.t()*solve(tV, tM);
	INTC(7, 0) = TEMP(0, 0);
	pr = TriNProb(wrap(tM), wrap(tV));
	if (pr < 0 || isnan(pr)) { pr = 0.0; } // rare numerical error
	TEMP = log(pr*sqrt(abs(det(tV)))*lam[0] * lam[0] * lam[0]) + 1.5*log(2.0*M_PI);
	INTC(7, 1) = TEMP(0, 0);

	// figure out the log Pbility  of the move
	arma::mat TT = INTC.col(1) - INTC.col(0);
	//cout << TT << endl; 
	arma::mat MT = max(TT);
	arma::mat SMULT = PMAT % exp(TT - repmat(MT, 8, 1));

	arma::mat log_prob = -0.5*(tau)*(Y.t()*Y) + log(sum(SMULT)) + MT;

	arma::mat POST_PROBS = exp(log(SMULT) - log(repmat(sum(SMULT), 8, 1))); //posterior sampling prob
	GetRNGstate();
	// get the cumulative  density
	double rvalue = R::runif(0, 1); int rval = -1;

	for (int i = 1; i < 8; i++) {
		POST_PROBS(i, 0) += POST_PROBS(i - 1, 0);

	}

	for (int i = 0; i < 8 && rval == -1; i++) {
		if (rvalue <= POST_PROBS(i, 0)) {
			rval = i;
		}
	}

	rval = rval + 1;

	PutRNGstate();
	arma::mat rVar, rMeans;

	switch (rval) {
	case 8:

		rVar = tV;
		rMeans = tM;
		break;
	case 1:

		rVar = 0.0;
		rMeans = 0.0;
		break;
	case 2:
	case 3:
	case 4:


		rVar = uniVars(rval - 2, 0);
		rMeans = uniMeans(rval - 2, 0);
		break;
	default:
		rVar = bivVars.submat(0, (rval - 5) * 2, 1, (rval - 5) * 2 + 1);
		rMeans = bivMeans.submat(0, rval - 5, 1, rval - 5);
		break;
	}
	List returnV(4);
	returnV[0] = wrap(rVar);
	returnV[1] = wrap(rMeans);
	returnV[2] = wrap(log_prob(0, 0));
	returnV[3] = wrap(rval);

	StringVector names(4);

	names(0) = "VAR";
	names(1) = "MEAN";
	names(2) = "LPROB";
	names(3) = "S_ELM";

	returnV.attr("names") = names;

	return returnV;
}

///////////////////////////////////////////////////////////////////////////
// sinsrtBeta_bivaraite and sdeleteBeta_bivariate are the same functions
// as the trivariate version just a different add or delete algorithm
//////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List sinsertBeta_bivariate(NumericVector tY, NumericVector Xidx, NumericVector ctau,
							NumericVector tp, NumericVector lam, NumericVector omega) {

	// initialize the proper values
	arma::mat Y = Rcpp::as<arma::mat>(tY);
	arma::mat X = Rcpp::as<arma::mat>(Xidx);

	arma::mat W = Rcpp::as<arma::mat>(omega);

	double tau = ctau[0];

	arma::mat tLAM = Rcpp::as<arma::mat>(lam);
	arma::mat LAM(tLAM.n_rows, tLAM.n_rows);  LAM.diag() = tLAM;

	arma::mat INTC(8, 2);

	INTC.zeros();
	arma::mat XW = X;
	XW.each_col() %= W;
	arma::mat B = (tau)*(XW.t()*X);          arma::mat U = (tau)*(XW.t()*Y) - lam[0];
	
	double p = tp[0];

	mat tV;            mat tM;
	mat uniMeans(3, 1); mat uniVars(3, 1);
	mat bivMeans(2, 3); mat bivVars(2, 6);

	mat TEMP; double pr = 0.0;


	mat PMAT(8, 1); PMAT(0, 0) = (p)*(p)*(p);

	// start with the univariate values

	for (int i = 1; i <4; i++) {
		PMAT(i, 0) = p * p*(1 - p);
		tV = (B(i - 1, i - 1));
		tV = tV.i();
		tM = tV * (U(i - 1, 0) - lam[0]);
		TEMP = -.5*tM.t()*solve(tV, tM);
		INTC(i, 0) = TEMP(0, 0);
		//Total positive area for the univariate normal distribution
		pr = 1 - R::pnorm(0, tM(0, 0), sqrt(tV(0, 0)), true, false);
		TEMP = log(pr*sqrt(tV)*lam[0]) + 0.5*log(2.0*M_PI);
		INTC(i, 1) = TEMP(0, 0);
		uniMeans(i - 1, 0) = tM(0, 0); uniVars(i - 1, 0) = tV(0, 0);

	}

	PMAT(4, 0) = p * (1 - p)*(1 - p);
	arma::mat tU = U.submat(0, 0, 1, 0);
	arma::mat tB = B.submat(0, 0, 1, 1);
	tV = (tB);
	tV = tV.i();
	tM = tV * (tU - lam[0]);
	bivVars.submat(0, 0, 1, 1) = tV;
	bivMeans.submat(0, 0, 1, 0) = tM;
	TEMP = -.5*tM.t()*solve(tV, tM);
	INTC(4, 0) = TEMP(0, 0);
	pr = BivNProb(wrap(tM), wrap(tV));
	if (pr < 0 || isnan(pr)) { pr = 0.0; } // rare numerical error
	TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0] * lam[0])) + log(2.0*M_PI);
	INTC(4, 1) = TEMP(0, 0);


	// Second Condition
	PMAT(5, 0) = p * (1 - p)*(1 - p);
	tB(0, 0) = B(0, 0); tB(0, 1) = B(0, 2);
	tB(1, 0) = B(2, 0); tB(1, 1) = B(2, 2);
	tU(0, 0) = U(0, 0); tU(1, 0) = U(2, 0);
	tV = (tB);
	tV = tV.i();
	tM = tV * (tU - lam[0]);
	bivVars.submat(0, 2, 1, 3) = tV;
	bivMeans.submat(0, 1, 1, 1) = tM;
	TEMP = -.5*tM.t()*solve(tV, tM);
	INTC(5, 0) = TEMP(0, 0);

	pr = BivNProb(wrap(tM), wrap(tV));
	if (pr < 0 || isnan(pr)) { pr = 0.0; } // rare numerical error
	TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0] * lam[0])) + log(2.0*M_PI);
	INTC(5, 1) = TEMP(0, 0);

	// Third Condition
	PMAT(6, 0) = p * (1 - p)*(1 - p);
	tU = U.submat(1, 0, 2, 0);
	tB = B.submat(1, 1, 2, 2);
	tV = (tB);
	tV = tV.i();
	tM = tV * (tU - lam[0]);
	bivVars.submat(0, 4, 1, 5) = tV;
	bivMeans.submat(0, 2, 1, 2) = tM;
	TEMP = -.5*tM.t()*solve(tV, tM);
	INTC(6, 0) = TEMP(0, 0);
	pr = BivNProb(wrap(tM), wrap(tV));
	if (pr < 0 || isnan(pr)) { pr = 0.0; } // rare numerical error
	TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0] * lam[0])) + log(2.0*M_PI);
	INTC(6, 1) = TEMP(0, 0);

	// TRIVARIATE CASE
	PMAT(7, 0) = (1 - p)*(1 - p)*(1 - p);
	tV = (B);
	tV = tV.i();
	tM = tV * (U - lam[0]);
	TEMP = -.5*tM.t()*solve(tV, tM);
	INTC(7, 0) = TEMP(0, 0);
	pr = TriNProb(wrap(tM), wrap(tV));
	if (pr < 0 || isnan(pr)) { pr = 0.0; } // rare numerical error
	TEMP = log(pr*sqrt(abs(det(tV)))*lam[0] * lam[0] * lam[0]) + 1.5*log(2.0*M_PI);
	INTC(7, 1) = TEMP(0, 0);

	// figure out the log Probility  of the move
	mat TT = INTC.col(1) - INTC.col(0);
	mat MT = max(TT);
	mat SMULT = PMAT % exp(TT - repmat(MT, 8, 1));
	mat log_prob = -0.5*(tau)*(Y.t()*Y) + log(sum(SMULT)) + MT;
	mat POST_PROBS = exp(log(SMULT) - log(repmat(sum(SMULT), 8, 1))); //posterior sampling prob
	GetRNGstate();

	// get the cumulative  density
	double rvalue = R::runif(0, 1); int rval = -1;

	for (int i = 1; i < 8; i++) {
		POST_PROBS(i, 0) += POST_PROBS(i - 1, 0);
	}

	for (int i = 0; i < 8 && rval == -1; i++) {
		if (rvalue <= POST_PROBS(i, 0)) {
			rval = i;
		}
	}

	rval = rval + 1;
	PutRNGstate();
	mat rVar, rMeans;

	int matDim = 2;
	switch (rval) {
	case 8:
		matDim = 3;
		rVar = tV;
		rMeans = tM;
		break;
	case 1:
		matDim = 1;
		rVar = 0.0;
		rMeans = 0.0;
		break;
	case 2:
	case 3:
	case 4:

		matDim = 1;
		rVar = uniVars(rval - 2, 0);
		rMeans = uniMeans(rval - 2, 0);
		break;
	default:

		rVar = bivVars.submat(0, (rval - 5) * 2, 1, (rval - 5) * 2 + 1);
		rMeans = bivMeans.submat(0, rval - 5, 1, rval - 5);
		break;
	}

	// COLLECT THE RESULTS AND PUT IT INTO A LIST  TO RETURN TO R
	List returnV(4);
	returnV[0] = wrap(rVar);
	returnV[1] = wrap(rMeans);
	returnV[2] = wrap(log_prob(0, 0));
	returnV[3] = wrap(rval);

	StringVector names(4);

	names(0) = "VAR";
	names(1) = "MEAN";
	names(2) = "LPROB";
	names(3) = "S_ELM";

	returnV.attr("names") = names;
	return returnV;
}


// [[Rcpp::export]]
List sdeleteBeta_bivariate(NumericVector tY, NumericVector Xidx, NumericVector ctau,
						   NumericVector tp, NumericVector lam, NumericVector omega) {

	// initialize the proper values
	arma::mat Y = Rcpp::as<arma::mat>(tY);
	arma::mat X = Rcpp::as<arma::mat>(Xidx);
	arma::mat W = Rcpp::as<arma::mat>(omega);
	double tau = ctau[0];

	arma::mat tLAM = Rcpp::as<arma::mat>(lam);
	arma::mat LAM(tLAM.n_rows, tLAM.n_rows);  LAM.diag() = tLAM;
	
	mat INTC(4, 2);

	INTC.zeros();
	arma::mat XW = X;
	 XW.each_col() %= W;
	arma::mat B = (tau)*(XW.t()*X);          arma::mat U = (tau)*(XW.t()*Y) - lam[0];
	double p = tp[0];

	mat tV;            mat tM;
	mat uniMeans(2, 1); mat uniVars(2, 1);
	mat TEMP; double pr = 0.0;


	mat PMAT(4, 1); PMAT(0, 0) = p * p;

	// start with the univariate values

	// start with the univariate values

	for (int i = 1; i <3; i++) {
		PMAT(i, 0) = p * (1 - p);
		tV = (B(i - 1, i - 1));
		tV = tV.i();
		tM = tV * (U(i - 1, 0) - lam[0]);
		TEMP = -.5*tM.t()*solve(tV, tM);
		INTC(i, 0) = TEMP(0, 0);
		//Total positive area for the univariate normal distribution
		pr = 1 - R::pnorm(0, tM(0, 0), sqrt(tV(0, 0)), true, false);
		TEMP = log(pr*sqrt(tV)*lam[0]) + 0.5*log(2.0*M_PI);
		INTC(i, 1) = TEMP(0, 0);
		uniMeans(i - 1, 0) = tM(0, 0); uniVars(i - 1, 0) = tV(0, 0);

	}

	//  bivariate values don't worry about looping as it just adds extra
	//  First  Condition

	PMAT(3, 0) = (1 - p)*(1 - p);
	mat tU = U;
	mat tB = B;
	tLAM = LAM;
	tV = (tB);
	tV = tV.i();
	tM = tV * (tU - lam[0]);

	TEMP = -.5*tM.t()*tV.i()*(tM);
	INTC(3, 0) = TEMP(0, 0);
	BivNProb(tM.memptr(), tV.memptr(), &pr);
	if (pr < 0) { pr = 0.0; } // rare numerical error
	TEMP = log(pr*(sqrt(abs(det(tV)))*sqrt(det(tLAM)))*4.0);

	INTC(3, 1) = TEMP(0, 0);

	// figure out the log Pbility  of the move

	mat TT = INTC.col(1) - INTC.col(0);
	mat MT = max(TT);
	mat SMULT = PMAT % exp(TT - repmat(MT, 4, 1));
	mat log_prob = -0.5*(tau)*(Y.t()*Y) + log(sum(SMULT)) + MT;
	mat POST_PROBS = exp(log(SMULT) - log(repmat(sum(SMULT), 4, 1))); //posterior sampling prob



	GetRNGstate();
	// get the cumulative  density
	double rvalue = R::runif(0, 1); int rval = -1;

	for (int i = 1; i < 4; i++) {
		POST_PROBS(i, 0) += POST_PROBS(i - 1, 0);
	}


	for (int i = 0; i < 4 && rval == -1; i++) {
		if (rvalue <= POST_PROBS(i, 0)) {
			rval = i;
		}
	}
	rval = rval + 1;
	PutRNGstate();


	mat rVar, rMeans;
	int matDim; 
	switch (rval) {
	case 4:
		matDim = 2;
		rVar = tV;
		rMeans = tM;
		break;
	case 1:
		matDim = 1;
		rVar = 0.0;
		rMeans = 0.0;
		break;
	default:
		matDim = 1;
		rVar = uniVars(rval - 2, 0);
		rMeans = uniMeans(rval - 2, 0);
		break;
	}

	// COLLECT THE RESULTS AND PUT IT INTO A LIST  TO RETURN TO R
	List returnV(4);
	returnV[0] = wrap(rVar);
	returnV[1] = wrap(rMeans);
	returnV[2] = wrap(log_prob(0, 0));
	returnV[3] = wrap(rval);

	StringVector names(4);

	names(0) = "VAR";
	names(1) = "MEAN";
	names(2) = "LPROB";
	names(3) = "S_ELM";

	returnV.attr("names") = names;
	return returnV;

}
