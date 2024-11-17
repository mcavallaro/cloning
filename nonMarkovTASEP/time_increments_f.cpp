#include <cstdio>
#include <vector>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "type_defs.h"
#include "functions.h"

double time_increment_f(System *system, int i, Parameters br, gsl_rng* r, vector<double> &times){
	int j;
    double beta = 0;
    double RV;
	int dim_chain = system[i].node.size();
	int N=0;
	for (j=0;j<dim_chain-1;j++){
		if (system[i].node[j] == 1){
			if (system[i].node[j+1] == 0){
				N++;
			}
		}
	}
	if (system[i].node[dim_chain-1]==1){
		beta = br.beta;
	}

	if (system[i].node[0]==1){
		return gsl_ran_exponential(r, 1/(N + beta ) ) ;
	}
	else if ((system[i].node[0]==0)&&(system[i].current==0)){
		return gsl_ran_exponential(r, 1/(N + beta + br.alpha0) ) ;
	}
	else {
//        printf("%f %f \n",times[0], times[i]);
		RV = ran_lambert(r,&system[0], i, br, N, beta, times);
		return RV;
	}
}
