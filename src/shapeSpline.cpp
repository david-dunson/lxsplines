#include <RcppArmadillo.h>
#include <iostream>

using namespace std;
using namespace arma; 
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector shapesplineInsertQuadratic(NumericVector k, NumericVector tp, NumericVector txi,
										 NumericVector tdeg, NumericVector tCBX, NumericVector tpos) {

	arma::mat knots = Rcpp::as<arma::mat>(k);
	arma::mat t = Rcpp::as<arma::mat>(tp);


	int deg = tdeg[0];
	int pos = tpos[0] - 1; // adjust  the position for C++ and extra knots
	knots.insert_rows(0, 1); knots.insert_rows(0, 1);
	knots.insert_rows(knots.n_rows, 1); knots.insert_rows(knots.n_rows, 1);
	knots(0, 0) = knots(2, 0); knots(1, 0) = knots(2, 0);
	knots(knots.n_rows - 1, 0) = knots(knots.n_rows - 3, 0);
	knots(knots.n_rows - 2, 0) = knots(knots.n_rows - 3, 0);

	IntegerVector arrayDims = tCBX.attr("dim");
	arma::cube CBX(tCBX.begin(), arrayDims[0], arrayDims[1], arrayDims[2]);
	arma::cube RVAL = arma::cube(arrayDims[0], arrayDims[1] + 1, arrayDims[2]);

	RVAL.zeros();
	//      Rprintf("%d %d %d %d",pos,TV.n_rows,TV.n_cols,TV.n_slices);
	RVAL.tube(0, 0, RVAL.n_rows - 1, pos - 1) = CBX.tube(0, 0, CBX.n_rows - 1, pos - 1);
	RVAL.tube(0, pos + 1, RVAL.n_rows - 1, RVAL.n_cols - 1) = CBX.tube(0, pos, CBX.n_rows - 1, CBX.n_cols - 1);

	for (int j = 0; j <= deg; j++) {
		for (int i = pos - 1; i <= pos + 2; i++) {
			double h1 = knots(i, 0); double h2 = knots(i + 1, 0);
			double h3 = knots(i + 2, 0); double h4 = knots(i + 3, 0);

			double    t1 = (((h2 - h1)*(h3 - h1)));
			double    t2 = (((h3 - h1)*(h3 - h2)));
			double    t3 = (((h4 - h2)*(h3 - h2)));
			double    t4 = (((h4 - h3)*(h4 - h2)));
			double    z1, z2, z3, z4;
			if (t1 == 0) { z1 = 0; }
			else { z1 = 1 / t1; }
			if (t2 == 0) { z2 = 0; }
			else { z2 = 1 / t2; }
			if (t3 == 0) { z3 = 0; }
			else { z3 = 1 / t3; }
			if (t4 == 0) { z4 = 0; }
			else { z4 = 1 / t4; }


			arma::umat LTH1 = (t <= h1);
			arma::umat LTH2 = (t <= h2);
			arma::umat LTH3 = (t <= h3);
			arma::umat LTH4 = (t <= h4);


			arma::mat S1 = 1.0 / double(3 + j)*pow(t, 3 + j) - 2.0*h1 / (2 + j)*pow(t, 2 + j) + pow(t, 1 + j)*pow(h1, 2) / double(1 + j); S1 = z1 * S1;
			arma::mat BS1 = S1 * 0 + (1.0 / double(3 + j)*pow(h1, 3 + j) - 2.0*h1 / double(2 + j)*pow(h1, 2 + j) + (pow(h1, 1 + j))*pow(h1, 2) / double(1 + j))*z1;
			arma::mat AS1 = S1 * 0 + (1.0 / double(3 + j)*pow(h2, 3 + j) - 2.0*h1 / double(2 + j)*pow(h2, 2 + j) + (pow(h2, 1 + j))*pow(h1, 2) / double(1 + j))*z1;

			arma::mat S2 = h3 / double(2 + j)*pow(t, 2 + j) - 1.0 / double(3 + j)*pow(t, 3 + j) - 1.0 / double(1 + j)*pow(t, j + 1)*h1*h3 + h1 / double(2 + j)*pow(t, 2 + j); S2 = S2 * z2;
			arma::mat BS2 = S2 * 0 + (h3 / double(2 + j)*pow(h2, 2 + j) - 1.0 / double(3 + j)*pow(h2, 3 + j) - 1.0 / double(1 + j)*pow(h2, j + 1)*h1*h3 + h1 / double(2 + j)*pow(h2, 2 + j))*z2;
			arma::mat AS2 = S2 * 0 + (h3 / double(2 + j)*pow(h3, 2 + j) - 1.0 / double(3 + j)*pow(h3, 3 + j) - 1.0 / double(1 + j)*pow(h3, j + 1)*h1*h3 + h1 / double(2 + j)*pow(h3, 2 + j))*z2;

			arma::mat S3 = h4 / double(2 + j)*pow(t, 2 + j) - 1.0 / double(3 + j)*pow(t, 3 + j) - h2 * h4 / double(1 + j)*pow(t, 1 + j) + h2 / double(2 + j)*pow(t, 2 + j); S3 = S3 * z3;
			arma::mat BS3 = S3 * 0 + (h4 / double(2 + j)*pow(h2, 2 + j) - 1.0 / double(3 + j)*pow(h2, 3 + j) - h2 * h4 / double(1 + j)*pow(h2, 1 + j) + h2 / double(2 + j)*pow(h2, 2 + j))*z3;
			AS2 = AS2 + (h4 / double(2 + j)*pow(h3, 2 + j) - 1.0 / double(3 + j)*pow(h3, 3 + j) - h2 * h4 / double(1 + j)*pow(h3, 1 + j) + h2 / double(2 + j)*pow(h3, 2 + j))*z3;

			arma::mat S4 = 1.0 / double(1 + j)*pow(h4, 2)*pow(t, 1 + j) - 2.0*h4 / double(2 + j)*pow(t, 2 + j) + 1.0 / double(3 + j)*pow(t, 3 + j); S4 = S4 * z4;
			arma::mat BS4 = S4 * 0 + (1.0 / double(1 + j)*pow(h4, 2)*pow(h3, 1 + j) - 2.0*h4 / double(2 + j)*pow(h3, 2 + j) + 1.0 / double(3 + j)*pow(h3, 3 + j))*z4;
			arma::mat AS3 = S4 * 0 + (1.0 / double(1 + j)*pow(h4, 2)*pow(h4, 1 + j) - 2.0*h4 / double(2 + j)*pow(h4, 2 + j) + 1.0 / double(3 + j)*pow(h4, 3 + j))*z4;

			arma::mat A = S1 * 0;
			A = A + ((S1 - BS1) % LTH2 + (AS1 - BS1) % (1 - LTH2)) % (1 - LTH1);
			A = A + ((S2 + S3 - (BS2 + BS3)) % LTH3 + (AS2 - (BS2 + BS3)) % (1 - LTH3)) % (1 - LTH2);
			A = A + ((S4 - BS4) % LTH4 + (AS3 - BS4) % (1 - LTH4)) % (1 - LTH3);
			RVAL.slice(j).submat(arma::span(), arma::span(i)) = A;
		}
	}


	return wrap(RVAL);

	// return R_NilValue;

}

// [[Rcpp::export]]
NumericVector shapesplineDeleteQuadratic(NumericVector k, NumericVector tp, NumericVector txi,
						NumericVector tdeg, NumericVector tCBX, NumericVector tpos) {
	arma::mat knots = Rcpp::as<arma::mat>(k);
	arma::mat t = Rcpp::as<arma::mat>(tp);

	int deg = tdeg[0];
	int pos = tpos[0] - 1; // adjust  the position for C++ indexing

	knots.insert_rows(0, 1); knots.insert_rows(0, 1);
	knots.insert_rows(knots.n_rows, 1); knots.insert_rows(knots.n_rows, 1); 
	knots(0, 0) = knots(2, 0); knots(1, 0) = knots(2, 0);
	knots(knots.n_rows - 1, 0) = knots(knots.n_rows - 3, 0);
	knots(knots.n_rows - 2, 0) = knots(knots.n_rows - 3, 0);

	IntegerVector arrayDims = tCBX.attr("dim");
	arma::cube CBX(tCBX.begin(), arrayDims[0], arrayDims[1], arrayDims[2]);
	arma::cube RVAL = arma::cube(arrayDims[0], arrayDims[1] - 1, arrayDims[2]);

	RVAL.zeros();
	//      Rprintf("%d %d %d %d",pos,TV.n_rows,TV.n_cols,TV.n_slices);
	RVAL.tube(0, 0, RVAL.n_rows - 1, pos - 1) = CBX.tube(0, 0, CBX.n_rows - 1, pos - 1);
	RVAL.tube(0, pos, RVAL.n_rows - 1, RVAL.n_cols - 1) = CBX.tube(0, pos + 1, CBX.n_rows - 1, CBX.n_cols - 1);

	for (int j = 0; j <= deg; j++) {
		for (int i = pos - 1; i <= pos + 1; i++) {
			double h1 = knots(i, 0); double h2 = knots(i + 1, 0);
			double h3 = knots(i + 2, 0); double h4 = knots(i + 3, 0);

			double    t1 = (((h2 - h1)*(h3 - h1)));
			double    t2 = (((h3 - h1)*(h3 - h2)));
			double    t3 = (((h4 - h2)*(h3 - h2)));
			double    t4 = (((h4 - h3)*(h4 - h2)));
			double    z1, z2, z3, z4;
			if (t1 == 0) { z1 = 0; }
			else { z1 = 1 / t1; }
			if (t2 == 0) { z2 = 0; }
			else { z2 = 1 / t2; }
			if (t3 == 0) { z3 = 0; }
			else { z3 = 1 / t3; }
			if (t4 == 0) { z4 = 0; }
			else { z4 = 1 / t4; }


			arma::umat LTH1 = (t <= h1);
			arma::umat LTH2 = (t <= h2);
			arma::umat LTH3 = (t <= h3);
			arma::umat LTH4 = (t <= h4);


			arma::mat S1 = 1.0 / double(3 + j)*pow(t, 3 + j) - 2.0*h1 / (2 + j)*pow(t, 2 + j) + pow(t, 1 + j)*pow(h1, 2) / double(1 + j); S1 = z1 * S1;
			arma::mat BS1 = S1 * 0 + (1.0 / double(3 + j)*pow(h1, 3 + j) - 2.0*h1 / double(2 + j)*pow(h1, 2 + j) + (pow(h1, 1 + j))*pow(h1, 2) / double(1 + j))*z1;
			arma::mat AS1 = S1 * 0 + (1.0 / double(3 + j)*pow(h2, 3 + j) - 2.0*h1 / double(2 + j)*pow(h2, 2 + j) + (pow(h2, 1 + j))*pow(h1, 2) / double(1 + j))*z1;

			arma::mat S2 = h3 / double(2 + j)*pow(t, 2 + j) - 1.0 / double(3 + j)*pow(t, 3 + j) - 1.0 / double(1 + j)*pow(t, j + 1)*h1*h3 + h1 / double(2 + j)*pow(t, 2 + j); S2 = S2 * z2;
			arma::mat BS2 = S2 * 0 + (h3 / double(2 + j)*pow(h2, 2 + j) - 1.0 / double(3 + j)*pow(h2, 3 + j) - 1.0 / double(1 + j)*pow(h2, j + 1)*h1*h3 + h1 / double(2 + j)*pow(h2, 2 + j))*z2;
			arma::mat AS2 = S2 * 0 + (h3 / double(2 + j)*pow(h3, 2 + j) - 1.0 / double(3 + j)*pow(h3, 3 + j) - 1.0 / double(1 + j)*pow(h3, j + 1)*h1*h3 + h1 / double(2 + j)*pow(h3, 2 + j))*z2;

			arma::mat S3 = h4 / double(2 + j)*pow(t, 2 + j) - 1.0 / double(3 + j)*pow(t, 3 + j) - h2 * h4 / double(1 + j)*pow(t, 1 + j) + h2 / double(2 + j)*pow(t, 2 + j); S3 = S3 * z3;
			arma::mat BS3 = S3 * 0 + (h4 / double(2 + j)*pow(h2, 2 + j) - 1.0 / double(3 + j)*pow(h2, 3 + j) - h2 * h4 / double(1 + j)*pow(h2, 1 + j) + h2 / double(2 + j)*pow(h2, 2 + j))*z3;
			AS2 = AS2 + (h4 / double(2 + j)*pow(h3, 2 + j) - 1.0 / double(3 + j)*pow(h3, 3 + j) - h2 * h4 / double(1 + j)*pow(h3, 1 + j) + h2 / double(2 + j)*pow(h3, 2 + j))*z3;

			arma::mat S4 = 1.0 / double(1 + j)*pow(h4, 2)*pow(t, 1 + j) - 2.0*h4 / double(2 + j)*pow(t, 2 + j) + 1.0 / double(3 + j)*pow(t, 3 + j); S4 = S4 * z4;
			arma::mat BS4 = S4 * 0 + (1.0 / double(1 + j)*pow(h4, 2)*pow(h3, 1 + j) - 2.0*h4 / double(2 + j)*pow(h3, 2 + j) + 1.0 / double(3 + j)*pow(h3, 3 + j))*z4;
			arma::mat AS3 = S4 * 0 + (1.0 / double(1 + j)*pow(h4, 2)*pow(h4, 1 + j) - 2.0*h4 / double(2 + j)*pow(h4, 2 + j) + 1.0 / double(3 + j)*pow(h4, 3 + j))*z4;

			arma::mat A = S1 * 0;
			A = A + ((S1 - BS1) % LTH2 + (AS1 - BS1) % (1 - LTH2)) % (1 - LTH1);
			A = A + ((S2 + S3 - (BS2 + BS3)) % LTH3 + (AS2 - (BS2 + BS3)) % (1 - LTH3)) % (1 - LTH2);
			A = A + ((S4 - BS4) % LTH4 + (AS3 - BS4) % (1 - LTH4)) % (1 - LTH3);
			RVAL.slice(j).submat(arma::span(), arma::span(i)) = A;
		}
	}

	return wrap(RVAL);

}

// [[Rcpp::export]]
NumericVector shapesplineInsertLinear(NumericVector k, NumericVector tp,
	NumericVector tdeg, NumericVector tCBX, NumericVector tpos) {

	arma::mat knots = Rcpp::as<arma::mat>(k);
	arma::mat t = Rcpp::as<arma::mat>(tp);



	int deg = tdeg[0];
	int pos = tpos[0] - 1; // adjust  the position for C++ indexing
	knots.insert_rows(0, 1); knots.insert_rows(knots.n_rows, 1);
	knots(0, 0) = knots(1, 0);
	knots(knots.n_rows - 1, 0) = knots(knots.n_rows - 2, 0);


	IntegerVector dims = tCBX.attr("dim");
	cube CBX = cube(tCBX.begin(), dims[0], dims[1], dims[2]);
	cube RVAL = cube(dims[0], dims[1] + 1, dims[2]);

	RVAL.zeros();
	// Rprintf("%d %d %d %d",pos,TV.n_rows,TV.n_cols,TV.n_slices);
	RVAL.tube(0, 0, RVAL.n_rows - 1, pos - 1) = CBX.tube(0, 0, CBX.n_rows - 1, pos - 1);
	RVAL.tube(0, pos + 1, RVAL.n_rows - 1, RVAL.n_cols - 1) = CBX.tube(0, pos, CBX.n_rows - 1, CBX.n_cols - 1);

	//cout << "The Knots :" << knots << endl; 

	for (int j = 0; j <= deg; j++) {

		for (int i = pos - 1; i <= pos + 1; i++) {
			double h1 = knots(i, 0); double h2 = knots(i + 1, 0); double h3 = knots(i + 2, 0);

			double t1 = ((h3 - h1)*(h2 - h1));
			double t2 = ((h3 - h1)*(h3 - h2));

			double z1 = (t1 == 0.0 ? 0.0 : 2.0 / t1);
			double z2 = (t2 == 0.0 ? 0.0 : 2.0 / t2);


			umat LTH1 = (t <= h1);
			umat LTH2 = (t <= h2);
			umat LTH3 = (t <= h3);


			double S1 = z1 * (1.0 / double(j + 2)* pow(h2, j + 2)
				- (h1 / double(j + 1))*pow(h2, j + 1)
				+ (1.0 / double(j + 1)
					- 1.0 / double(j + 2))*pow(h1, j + 2));
			double S2 = z2 * (h3 / double(j + 1)*pow(h3, j + 1)
				- 1.0 / double(j + 2)*pow(h3, (j + 2))
				+ (pow(h2, j + 2) / double(j + 2)
					- h3 * pow(h2, j + 1) / double(j + 1)));

			mat  A1 = z1 * (1.0 / double(j + 2)*pow(t, j + 2)
				- h1 / double(j + 1)*pow(t, j + 1)
				+ (1.0 / double(j + 1)
					- 1.0 / double(j + 2))*pow(h1, j + 2));
			mat  A2 = S1 + z2 * (h3 / double(j + 1)*pow(t, j + 1) - 1.0 / double(j + 2)*pow(t, j + 2) + (pow(h2, j + 2) / double(j + 2) - h3 * pow(h2, j + 1) / double(j + 1)));


			mat TEMP = conv_to<mat>::from((1 - LTH1) % (1 - LTH2) % (1 - LTH3));

			RVAL.slice(j).submat(span(), span(i)) = ((1 - LTH1) % LTH2%A1 + (1 - LTH1) % (1 - LTH2) % LTH3%A2 + TEMP * (S1 + S2));

		}
	}

	return wrap(RVAL);

}

// [[Rcpp::export]]
NumericVector shapesplineDeleteLinear(NumericVector k, NumericVector tp,
						NumericVector tdeg, NumericVector tCBX, NumericVector tpos) {
	arma::mat knots = Rcpp::as<arma::mat>(k);
	arma::mat t = Rcpp::as<arma::mat>(tp);

	int deg = tdeg[0];
	int pos = tpos[0] - 1; // adjust  the position for C++ indexing

	knots.insert_rows(0, 1); knots.insert_rows(knots.n_rows, 1);
	knots(0, 0) = knots(1, 0);
	knots(knots.n_rows - 1, 0) = knots(knots.n_rows - 2, 0);

	IntegerVector dims = tCBX.attr("dim");
	cube CBX = cube(tCBX.begin(), dims[0], dims[1], dims[2]);
	cube RVAL = cube(dims[0], dims[1] - 1, dims[2]);

	RVAL.zeros();
	//      Rprintf("%d %d %d %d",pos,TV.n_rows,TV.n_cols,TV.n_slices);
	RVAL.tube(0, 0, RVAL.n_rows - 1, pos - 1) = CBX.tube(0, 0, CBX.n_rows - 1, pos - 1);
	RVAL.tube(0, pos, RVAL.n_rows - 1, RVAL.n_cols - 1) = CBX.tube(0, pos + 1, CBX.n_rows - 1, CBX.n_cols - 1);



	for (int j = 0; j <= deg; j++) {
		for (int i = pos - 1; i <= pos; i++) {
			double h1 = knots(i, 0); double h2 = knots(i + 1, 0); double h3 = knots(i + 2, 0);

			double t1 = ((h3 - h1)*(h2 - h1));
			double t2 = ((h3 - h1)*(h3 - h2));

			double z1 = (t1 == 0.0 ? 0.0 : 2.0 / t1);
			double z2 = (t2 == 0.0 ? 0.0 : 2.0 / t2);


			umat LTH1 = (t <= h1);
			umat LTH2 = (t <= h2);
			umat LTH3 = (t <= h3);


			double S1 = z1 * (1.0 / double(j + 2)* pow(h2, j + 2)
				- (h1 / double(j + 1))*pow(h2, j + 1)
				+ (1.0 / double(j + 1)
					- 1.0 / double(j + 2))*pow(h1, j + 2));
			double S2 = z2 * (h3 / double(j + 1)*pow(h3, j + 1)
				- 1.0 / double(j + 2)*pow(h3, (j + 2))
				+ (pow(h2, j + 2) / double(j + 2)
					- h3 * pow(h2, j + 1) / double(j + 1)));

			mat  A1 = z1 * (1.0 / double(j + 2)*pow(t, j + 2)
				- h1 / double(j + 1)*pow(t, j + 1)
				+ (1.0 / double(j + 1)
					- 1.0 / double(j + 2))*pow(h1, j + 2));
			mat  A2 = S1 + z2 * (h3 / double(j + 1)*pow(t, j + 1) - 1.0 / double(j + 2)*pow(t, j + 2) + (pow(h2, j + 2) / double(j + 2) - h3 * pow(h2, j + 1) / double(j + 1)));


			mat TEMP = conv_to<mat>::from((1 - LTH1) % (1 - LTH2) % (1 - LTH3));

			RVAL.slice(j).submat(span(), span(i)) = ((1 - LTH1) % LTH2%A1 + (1 - LTH1) % (1 - LTH2) % LTH3%A2 + TEMP * (S1 + S2));

		}
	}

	return wrap(RVAL);

}
