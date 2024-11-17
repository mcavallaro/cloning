#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include "type_defs.h"
#include "time_jumps.h"

int assign_part_sum(Systems *system,int i);
int update(System*, int, int );
inline int rand_to_int(double, gsl_rng*);
double time_increment_f(System *system, int i, gsl_rng* r, vector<double> &, vector<double> &);

void usage(int argc, char *argv[]){
  printf("***************************************************************\n"
         "**                                                           **\n"
         "**   -*-   SemiMarkov ASEP MODEL - ONE NODE CLONING   -*-    **\n"
         "**                                                           **\n"
         "***************************************************************\n"
         "\n\n"
         " This is Free Software - You can use and distribute it under \n"
         " the terms of the GNU General Public License, version 3 or later\n\n"
         " (c) Massimo Cavallaro (m.cavallaro@qmul.ac.uk)\n\n");
  printf("Usage: %s [# simulation time] [DIM. ENSEMBLE] [s] [folder]\nargc=%d\n" , argv[0], argc);
}

int reservoir_sampling(gsl_rng * r,std::vector<int> &choices, int N, int k ) {
    if (k > N) return EXIT_FAILURE;
    int j,i;
    choices.clear();
    for (i = 0; i < k; i++){
        choices.push_back(i);
    }
    for (; i < N; i++) {
        j = gsl_rng_uniform_int(r, i+1);
        if (j < k) {
          choices[j] = i;
        }
    }
    return EXIT_SUCCESS;
}


int main(int argc, char *argv[]) {

	if (argc != 5 ){
	    usage(argc,argv);
	    exit(1);
	}

	FILE* fout;

	int i;
	int DIM_ENSEMBLE;
	double T;
	double s = atof(argv[3]);
	int new_copies;
	int num_rand;
	int y;
	int j;

	double cloning;

	double R;
	double clfact;

	gsl_rng * r;
	double mean_observable;

	int status = 0;

	T = atof(argv[1]);
	DIM_ENSEMBLE = atoi(argv[2]);

	char nomefile [80];

	vector<int> choice;
	vector<int>::iterator choiceI;
	vector<double> jump_time(DIM_ENSEMBLE);

	double next_time;

    vector<double> CDF1;
    vector<double> CDF0;
    CDF0.resize(1000000);
    CDF1.resize(1000000);

    double bin_size = 0.001;

    double cdf = RTD_1(bin_size*0.5);
    int ssize = CDF1.size();
    CDF1[0] = cdf;
    for (j=1; (j<ssize )&&(cdf < 0.99999); j++){
        cdf = RTD_1(bin_size*(j+0.5));
        CDF1[j] =cdf;
    }
    CDF1[j]=1;
    CDF1.resize(j);

    cdf = RTD_0(bin_size*0.5);
    ssize = CDF0.size();
    CDF0[0] = cdf;
    for (j=1; (j<ssize) &&(cdf < 0.99999); j++){
        cdf = RTD_0(bin_size*(j+0.5));
        CDF0[j] =cdf;
    }
    CDF0[j]=1;
    CDF0.resize(j);



	vector<double> clfactv(4,1);
	clfactv[2]=exp(s);
	clfactv[3]=exp(-s);

	vector<Systems> system(DIM_ENSEMBLE);
	choice.reserve(10000000);

	r = gsl_rng_alloc(gsl_rng_taus);

	sprintf(nomefile,"%sASEP%s_%s_%s.txt", argv[4],argv[1],argv[2],argv[3] );
	fout = fopen(nomefile,"w");

	for (i=0; i<DIM_ENSEMBLE; i++) {
	    system[i].part_sum.resize(4);
	    system[i].part_sum.assign(4,0);
	}

	for (i=0; i<DIM_ENSEMBLE; i++) {
		system[i].current=0;
		system[i].site=0;
	}

	for (i=0; i<DIM_ENSEMBLE; i++) {
		jump_time[i] = system[i].time = time_increment_f(&system[0], i, r, CDF0, CDF1);;
	}

	cloning = 0;

    i = min_element(jump_time.begin(), jump_time.end()) - jump_time.begin();

	while (jump_time[i] < T){
        status = assign_part_sum(&system[0],i);
        if (status != EXIT_SUCCESS){
            fprintf(stderr,"PROBLEM in assign_part_sum %d", i);
	        return EXIT_FAILURE;
        }

	    R = gsl_rng_uniform(r) * system[i].part_sum.back();
	    j = upper_bound(system[i].part_sum.begin(),system[i].part_sum.end(),R) - system[i].part_sum.begin();
	    status = update(&system[0], i, j);
	    if ( status != 1 ) {
	    	cerr << "ERROR " << i << " " << j << endl;
	        return EXIT_FAILURE;
	    }

	    next_time = time_increment_f(&system[0], i, r, CDF0, CDF1);

	    jump_time[i] = jump_time[i] + next_time;
        system[i].time = next_time; 

       /* Compute the cloning factor */
	    clfact  = clfactv[j];
	    if (clfact != 1) {
	        y = rand_to_int(clfact,r);
	    }
	    else {
	        y=1;
	    }

	    if (y == 0) {
	    	// do{num_rand = gsl_rng_uniform_int(r,DIM_ENSEMBLE);}while(num_rand==i);
            num_rand = gsl_rng_uniform_int(r, DIM_ENSEMBLE - 1);
            if (num_rand >= i) num_rand++;
	        system[i] = system[num_rand];
	        jump_time[i] = jump_time[num_rand];
	    }
	    else if (y > 1) {
		    new_copies = y - 1;
		    status = reservoir_sampling(r,choice, DIM_ENSEMBLE+new_copies, new_copies);
		    if (status!=EXIT_SUCCESS) {
		    	cerr << "ERROR in reserve_sampling" << endl;
		        return EXIT_FAILURE;
		    }
		    for (choiceI = choice.begin(); choiceI != choice.end() ; choiceI++) {
			    if ((*choiceI) < DIM_ENSEMBLE) {
				    system[*choiceI] = system[i];
			    	jump_time[*choiceI]=jump_time[i];
			    }
			}
        }

	    cloning = cloning + log( (double)(DIM_ENSEMBLE+clfact-1)/(double)DIM_ENSEMBLE);	     
	    i = min_element(jump_time.begin(), jump_time.end()) - jump_time.begin();
	}

    fprintf(fout,"%f\t%f\t", s, cloning/T);

	/* Thermodynamic integration */
	mean_observable=0;
	for (int it=0; it<DIM_ENSEMBLE; it++ ) {
		mean_observable = mean_observable + (double)(system[it].current)/T;
	}
	fprintf(fout,"%e\n", mean_observable/(double)DIM_ENSEMBLE);



    double SUM = 0;
    for (int it=0; it<DIM_ENSEMBLE; it++){
        SUM = SUM + system[it].site;
    }
    printf("prob: %f %f\n", SUM/DIM_ENSEMBLE, (DIM_ENSEMBLE - SUM)/DIM_ENSEMBLE);

	gsl_rng_free(r);
	fclose(fout);
	return EXIT_SUCCESS;
}

inline int rand_to_int(double clfact, gsl_rng * r){
    return    (clfact + gsl_rng_uniform(r));
}

