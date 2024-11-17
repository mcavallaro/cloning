#include <cstdlib>
#include <gsl/gsl_rng.h>
#include "type_defs.h"



#ifndef TIME_JUMPS
#define TIME_JUMPS

double RTD_0(double tau);
double RTD_1(double tau);


double Palpha(Systems *system, int i);
double Pbeta(Systems *system, int i);
double Pgamma(Systems *system, int i);
double Pdelta(Systems *system, int i);

#endif 


