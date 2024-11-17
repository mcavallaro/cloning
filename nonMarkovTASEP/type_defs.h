#include <vector>

using namespace std;

#ifndef type_def
#define type_def
typedef struct Systems {
	vector<int> node;
	double current;	
	vector <double> part_sum;
	int N;
//	double total_time;
//	double history;
	double age;
} System;


typedef struct Param {
	double a;
	double alpha0;
    double beta;
} Parameters;
#endif
