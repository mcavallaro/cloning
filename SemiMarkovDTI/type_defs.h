#include <cstdlib>
#include <vector>

using namespace std;

typedef struct Systems {
	int site;
	int current;
	double time_increment;
	vector <double> part_sum;
/*	vector <double> part_sum_conserved;*/
} System;



typedef struct probabilities {
	double alpha;
	double beta;
	double gamma;
	double delta;
} Probabilities;






