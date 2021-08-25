# include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double llh_fun(NumericVector Y, NumericVector d, 
	NumericMatrix X0, NumericMatrix X1, NumericMatrix X2,
	NumericVector alpha, NumericVector beta1, NumericVector beta2, double t0) {

	int n = Y.size();
	int p0 = alpha.size();
	int p1 = beta1.size();
	int p2 = beta2.size();

	NumericVector theta0t(n);
	NumericVector theta1t(n);
	NumericVector theta2t(n);

	NumericVector theta0(n);
	NumericVector theta1(n);
	NumericVector theta2(n);

	NumericVector temp1(n);
	NumericVector temp2(n);
	NumericVector tempcheck(n);
	double out;

	for (int i = 0; i < n; i ++) {
		for (int j0 = 0; j0 < p0; j0 ++) {
			theta0t[i] += X0(i, j0) * alpha[j0];
		}

		theta0[i] = exp(theta0t[i]);

		for (int j1 = 0; j1 < p1; j1 ++) {
			theta1t[i] += X1(i, j1) * beta1[j1];
		}
		theta1[i] = exp(theta1t[i]);

		for (int j2 = 0; j2 < p2; j2 ++) {
			theta2t[i] += X2(i, j2) * beta2[j2];
		}		
		theta2[i] = exp(theta2t[i]);

		tempcheck[i] = 1 - (theta0[i] * Y[i] + 
			theta1[i] * (pow(Y[i]+t0, 1 - theta2[i]) / (1-theta2[i]) - 
				pow(t0, 1-theta2[i]) / (1-theta2[i])));
		//std::cout<<theta2[i];
		if (tempcheck[i] < 0 | tempcheck[i] > 1) {
			out = -999999999;
			return(out);
		}

		if (d[i] == 0) {
			temp2[i] = log(tempcheck[i]);
			temp1[i] = 0;
		} else {
			temp2[i] = 0;
			temp1[i] = log(theta0[i] + theta1[i] * pow(Y[i]+t0, -theta2[i]));
		}

	}

	out = sum(temp1) + sum(temp2);
	return out;
}



// [[Rcpp::export]]
double llh_fun_nomiss(NumericVector Y, NumericVector d, 
	NumericMatrix X0, NumericVector alpha) {

	int n = Y.size();
	int p0 = alpha.size();

	NumericVector theta0t(n);

	NumericVector theta0(n);

	NumericVector temp1(n);
	NumericVector temp2(n);
	NumericVector tempcheck(n);
	double out;

	for (int i = 0; i < n; i ++) {
		for (int j0 = 0; j0 < p0; j0 ++) {
			theta0t[i] += X0(i, j0) * alpha[j0];
		}

		theta0[i] = exp(theta0t[i]);

		tempcheck[i] = 1 - (theta0[i] * Y[i]);
		if (tempcheck[i] < 0 | tempcheck[i] > 1) {
			out = -999999999;
			return(out);
		}

		if (d[i] == 0) {
			temp2[i] = log(tempcheck[i]);
			temp1[i] = 0;
		} else {
			temp2[i] = 0;
			temp1[i] = log(theta0[i]);
		}

	}

	out = sum(temp1) + sum(temp2);
	return out;
}


// [[Rcpp::export]]
double llh_fun_exp(NumericVector Y, NumericVector d, 
	NumericMatrix X0, NumericMatrix X1, NumericMatrix X2,
	NumericVector alpha, NumericVector beta1, NumericVector beta2) {

	int n = Y.size();
	int p0 = alpha.size();
	int p1 = beta1.size();
	int p2 = beta2.size();

	NumericVector theta0t(n);
	NumericVector theta1t(n);
	NumericVector theta2t(n);

	NumericVector theta0(n);
	NumericVector theta1(n);
	NumericVector theta2(n);

	NumericVector temp1(n);
	NumericVector temp2(n);
	NumericVector tempcheck(n);
	double out;

	for (int i = 0; i < n; i ++) {
		for (int j0 = 0; j0 < p0; j0 ++) {
			theta0t[i] += X0(i, j0) * alpha[j0];
		}

		theta0[i] = exp(theta0t[i]);

		for (int j1 = 0; j1 < p1; j1 ++) {
			theta1t[i] += X1(i, j1) * beta1[j1];
		}
		theta1[i] = exp(theta1t[i]);

		for (int j2 = 0; j2 < p2; j2 ++) {
			theta2t[i] += X2(i, j2) * beta2[j2];
		}		
		theta2[i] = exp(theta2t[i]) / (1 + exp(theta2t[i]));

		tempcheck[i] = 1 - (theta0[i] * Y[i] + 
			theta1[i] / log(theta2[i]) * (pow(theta2[i], Y[i]) - 1));
		if (tempcheck[i] < 0 | tempcheck[i] > 1) {
			out = -999999999;
			return(out);
		}

		if (d[i] == 0) {
			temp2[i] = log(tempcheck[i]);
			temp1[i] = 0;
		} else {
			temp2[i] = 0;
			temp1[i] = log(theta0[i] + theta1[i] * pow(theta2[i], Y[i]));
		}

	}

	out = sum(temp1) + sum(temp2);
	return out;
}





// [[Rcpp::export]]
double llh_fun_poly_discrete(NumericVector Y, NumericVector d, 
	NumericMatrix X0, NumericMatrix X1, NumericMatrix X2,
	NumericVector alpha, NumericVector beta1, NumericVector beta2, double t0) {

	int n = Y.size();
	int p0 = alpha.size();
	int p1 = beta1.size();
	int p2 = beta2.size();

	NumericVector theta0t(n);
	NumericVector theta1t(n);
	NumericVector theta2t(n);

	NumericVector theta0(n);
	NumericVector theta1(n);
	NumericVector theta2(n);

	NumericVector temp1(n);
	NumericVector temp2(n);
	NumericVector tempcheck(n);
	double out;

	for (int i = 0; i < n; i ++) {
		for (int j0 = 0; j0 < p0; j0 ++) {
			theta0t[i] += X0(i, j0) * alpha[j0];
		}

		theta0[i] = exp(theta0t[i]);

		for (int j1 = 0; j1 < p1; j1 ++) {
			theta1t[i] += X1(i, j1) * beta1[j1];
		}
		theta1[i] = exp(theta1t[i]);

		for (int j2 = 0; j2 < p2; j2 ++) {
			theta2t[i] += X2(i, j2) * beta2[j2];
		}		
		theta2[i] = exp(theta2t[i]);

		tempcheck[i] = 1 - (theta0[i] * ((Y[i]-1>0?(Y[i]-1):0)) + 
			theta1[i] * (pow((Y[i]+t0-1>t0?(Y[i]+t0-1):t0), 1 - theta2[i]) / (1-theta2[i]) - 
				pow(t0, 1-theta2[i]) / (1-theta2[i])));
		if (tempcheck[i] < 0 | tempcheck[i] > 1) {
			out = -999999999;
			return(out);
		}

		if (d[i] == 0) {
			temp2[i] = log(tempcheck[i]);
			temp1[i] = 0;
		} else {
			temp2[i] = 0;
			temp1[i] = log((Y[i]>1?(theta0[i]*2):(theta0[i])) + theta1[i] * (pow(Y[i]+t0 + 1, 1-theta2[i]) - 
				pow((Y[i]>1?(Y[i]+t0 - 1):t0+1), 1 - theta2[i])) / (1-theta2[i]));
		}

	}

	out = sum(temp1) + sum(temp2);
	return out;
}


// [[Rcpp::export]]
double llh_fun_nomis_discrete(NumericVector Y, NumericVector d, 
	NumericMatrix X0, NumericVector alpha) {

	int n = Y.size();
	int p0 = alpha.size();

	NumericVector theta0t(n);	
	NumericVector theta0(n);


	NumericVector temp1(n);
	NumericVector temp2(n);
	NumericVector tempcheck(n);
	double out;

	for (int i = 0; i < n; i ++) {
		for (int j0 = 0; j0 < p0; j0 ++) {
			theta0t[i] += X0(i, j0) * alpha[j0];
		}

		theta0[i] = exp(theta0t[i]);

		tempcheck[i] = 1 - (theta0[i] * (Y[i]-1>0?(Y[i]-1):0));
		if (tempcheck[i] < 0 | tempcheck[i] > 1) {
			out = -999999999;
			return(out);
		}

		if (d[i] == 0) {
			temp2[i] = log(tempcheck[i]);
			temp1[i] = 0;
		} else {
			temp2[i] = 0;
			temp1[i] = log(theta0[i] * 2);
		}

	}

	out = sum(temp1) + sum(temp2);
	return out;
}



// [[Rcpp::export]]
double llh_fun_extpoly(NumericVector Y, NumericVector d, 
	NumericMatrix X0, NumericMatrix X1, NumericMatrix X2, NumericMatrix X3,
	NumericVector alpha, NumericVector beta1, NumericVector beta2, NumericVector beta3) {

	int n = Y.size();
	int p0 = alpha.size();
	int p1 = beta1.size();
	int p2 = beta2.size();
	int p3 = beta3.size();

	NumericVector theta0t(n);
	NumericVector theta1t(n);
	NumericVector theta2t(n);
	NumericVector theta3t(n);

	NumericVector theta0(n);
	NumericVector theta1(n);
	NumericVector theta2(n);
	NumericVector theta3(n);

	NumericVector temp1(n);
	NumericVector temp2(n);
	NumericVector tempcheck(n);
	double out;

	for (int i = 0; i < n; i ++) {
		for (int j0 = 0; j0 < p0; j0 ++) {
			theta0t[i] += X0(i, j0) * alpha[j0];
		}

		theta0[i] = exp(theta0t[i]);

		for (int j1 = 0; j1 < p1; j1 ++) {
			theta1t[i] += X1(i, j1) * beta1[j1];
		}
		theta1[i] = exp(theta1t[i]);

		for (int j2 = 0; j2 < p2; j2 ++) {
			theta2t[i] += X2(i, j2) * beta2[j2];
		}		
		theta2[i] = exp(theta2t[i]);

		for (int j3 = 0; j3 < p3; j3 ++) {
			theta3t[i] += X3(i, j3) * beta3[j3];
		}
		theta3[i] = exp(theta3t[i]);

		tempcheck[i] = 1 - (theta0[i] * Y[i] + 
			theta1[i] * (pow(Y[i]+theta3[i], 1 - theta2[i]) / (1-theta2[i]) - 
				pow(theta3[i], 1-theta2[i]) / (1-theta2[i])));
		if (tempcheck[i] < 0 | tempcheck[i] > 1) {
			out = -999999999;
			return(out);
		}

		if (d[i] == 0) {
			temp2[i] = log(tempcheck[i]);
			temp1[i] = 0;
		} else {
			temp2[i] = 0;
			temp1[i] = log(theta0[i] + theta1[i] * pow(Y[i]+theta3[i], -theta2[i]));
		}

	}

	out = sum(temp1) + sum(temp2);
	return out;
}
