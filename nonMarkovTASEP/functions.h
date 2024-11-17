#include <gsl/gsl_rng.h>

double w_0(System *system, int i, Parameters br, vector<double> &history);

double phi_c(Systems *system, int i,Parameters br, int N, double beta, double tau);

/*double w_j(System *system, int i, Parameters br, double timet);*/
double ran_lambert(gsl_rng* r, Systems *system, int i, Parameters br, int N, double beta, vector<double> &history);
