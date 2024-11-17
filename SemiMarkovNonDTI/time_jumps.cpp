#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "type_defs.h"
#include <gsl/gsl_sf_gamma.h>
#include <cmath>

//double gsl_ran_weibull_pdf (double x, double a, double b); //WTD
//double gsl_cdf_weibull_Q (double x, double a, double b); //survivor

// We first fix the Shape parameter
/*double k01L = 1.5;
double k10L = 1.5;
double k01R = 0.5;
double k10R = 0.5;*/

double k01L = 2;
double k10L = 2;
double k01R = 2;
double k10R = 2;

//In order to compare it to the Exponential case, we first fix the effective rates
double m01L = 0.1;
double m10L = 0.2;
double m01R = 0.3;
double m10R = 0.2;

//then we can derive the Scale parameters.
double alpha = m01L/k01L;
double Beta  = m10L/k10L;
double gamm  = m01R/k01R;
double delta = m10R/k10R;

/*
DO NOT MODIFY THE PART BELOW
*/
inline double hazard(double tau, double  lambda, double  k){
	return	exp(-tau/lambda) *pow(tau/lambda ,k) /tau /gsl_sf_gamma_inc(k,tau/lambda);
}

double RTD_0(double tau){ // duration of stay in state 0
	return 1 - gsl_cdf_gamma_Q(tau,k10L,Beta)*gsl_cdf_gamma_Q(tau,k10R,delta) ;
}

double RTD_1(double tau){ // duration of stay in state 1
	return 1 - gsl_cdf_gamma_Q(tau,k01L,alpha)*gsl_cdf_gamma_Q(tau,k01R,gamm) ;
}

double Palpha(Systems *system, int i){
	double tau = system[i].time;
    return  hazard(tau,alpha,k01L);	
}
double Pbeta(Systems *system, int i){
	double tau = system[i].time;
    return   hazard(tau,Beta,k10L);
}

double Pgamma(Systems *system, int i){
	double tau = system[i].time;
    return   hazard(tau,gamm,k01R);
}

double Pdelta(Systems *system, int i){
	double tau = system[i].time;
    return   hazard(tau,delta,k10R);
}

double time_increment_f( Systems* system,int i, gsl_rng* r, vector<double> &CDF0, vector<double> &CDF1){
	double R = gsl_rng_uniform(r);
	int g;

	if (system[i].site==1){
		g = upper_bound(CDF1.begin(),CDF1.end(), R)-CDF1.begin();		
/*		cout << system[i].site << " " << g << " " << (g+0.5)*0.00001 << endl;*/
	}
	else if (system[i].site==0){
		g = upper_bound(CDF0.begin(),CDF0.end(), R)-CDF0.begin();
/*		cout << system[i].site << " " << g << " " << (g+0.5)*0.00001 << endl;*/
	}
	else {
		fprintf(stderr, "%s\n", "minchia" );
		exit(EXIT_FAILURE);
	}
	// printf("%.2f\n", (g+0.5)*0.001);
	return (g+0.5)*0.001;
}