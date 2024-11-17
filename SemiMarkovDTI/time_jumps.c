#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "type_defs.h"


/*double k1 =1;  //SCGF1.txt
double theta1 =1 ;
double theta0 =1;
double k0 =1 ;*/

/*double theta0 = 100;  //SCGF2.txt
double k0 =1./10. ;
double k1 = 1 ;
double theta1 =1 ;*/


double theta0 = 100; //SCGF3.txt
double k0 =1. ;
double k1 = 1./10. ;
double theta1 =1 ;

double time_increment(System *system, int i, gsl_rng* r){
//    return gsl_ran_exponential(r,1./system[i].part_sum.back());
    if (system[i].site){
        return gsl_ran_gamma(r, k1, theta1);
    }
    else {
        return gsl_ran_gamma(r, k0, theta0);
    }
}
