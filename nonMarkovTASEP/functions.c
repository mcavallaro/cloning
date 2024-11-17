#include <cstdlib>
#include <iostream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include "type_defs.h"



/*double phi_c(Systems *system, int i, Parameters br, int N, double beta, double tau, vector<double> &history){
    double power = pow( history[i]/(history[i]+tau), br.a*system[i].current );
//    double expon = exp(-(  br.alpha0  + N + beta) * tau);
    double expon = exp(-(  br.alpha0 + N + beta) * tau);
    return 1 - power * expon;
}*/


inline double inv_deriv(double x){
	return -(x+1)/x;
//	return -1.-1./x;
}

double MyLambert(double exponent, double A){
	/*We want to find the solution of:
	 W0(A e^exponent)=y
	 A * e^exponent = y * e^y
	 ln (A) + exponent = ln(y) + y
	*/
	 double constant = log(A) + exponent;  // this is the first apporox for the solution y (as log(y)<<y ) 'constant' is also large than the solution 'y'.
	 double yold = constant;
	 double y;
	 double eps=99999;
//	 	printf(" %f %f %f \n",eps, exponent, A);
	 int n=0;
	 while((  eps > 0.0000001)&&(n<1000) ){
	 	y = yold - (constant - log(yold) - yold)*inv_deriv(yold) ;
	 	eps = fabs(yold - y);
	 	yold = y;
	 	n++;
//	 	printf(" %e \n",eps);
	 }
	 return y;
}



double ran_lambert(gsl_rng* r, Systems *system, int i, Parameters br, int N, double beta, vector<double> &history){

	double summ  = N + beta + br.alpha0;
	double R     = gsl_rng_uniform(r);
	double acurr = br.a*system[i].current;
//	printf("it %f %f %f %f\n", history[i], summ, R, system[i].current);
	double Wb = MyLambert(history[i] * summ/acurr, history[i]*summ*pow(1-R,-1/acurr) /acurr);
	if (history[i] > acurr*Wb/summ){
		printf("if  %f\n", -history[i] + acurr*Wb/summ);
	}
	return -history[i] + acurr*Wb/summ;
}


double w_0(Systems *system, int i, Parameters br, vector<double> &history){
    return br.alpha0  + br.a* system[i].current/history[i];
}





