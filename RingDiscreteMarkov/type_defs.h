#include <cstdlib>
#include <vector>

using namespace std;


#ifndef TYPE_DEF
#define TYPE_DEF

typedef struct Systems {
	int position;
	double current;
	vector<double> part_sum;
} System;

#endif
