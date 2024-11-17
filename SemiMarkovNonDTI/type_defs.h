#include <cstdlib>
#include <vector>

#ifndef TYPE_DEFS
#define TYPE_DEFS

using namespace std;

typedef struct Systems {
	int site;
	int current;
	double time;
	vector <double> part_sum;
} System;


#endif