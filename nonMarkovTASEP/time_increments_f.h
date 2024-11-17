#include <vector>
#include <gsl/gsl_rng.h>
#include "type_defs.h"


double time_increment_f(System *, int , Parameters , gsl_rng* , vector<double> &);
