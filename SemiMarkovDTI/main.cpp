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

int update(System*, int, int, Probabilities, double );
inline int rand_to_int(double, gsl_rng*);
double time_increment(System *system, int i, gsl_rng* r);

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
  printf("Usage: %s [# simulation time] [DIM. ENSEMBLE] [alpha] [beta] [gamma] [delta] [s] [folder]\nargc=%d\n" , argv[0], argc);
  printf("Warning -- [alpha] [beta] [gamma] and [delta] must sum to one.\n");
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

	if (argc != 9 ){
	    usage(argc,argv);
	    exit(1);
	}

	FILE* fout;

	int i;
	int T, DIM_ENSEMBLE;
	double s = atof(argv[7]);
	int new_copies;
	int num_rand;
	int y;
	int j;

	Probabilities Pr;

	long double cloning;

	double R;
	double clfact;

	gsl_rng * r;
	double mean_observable;

	int status = 0;
	int weakness, steps;

	T = atof(argv[1]);
	DIM_ENSEMBLE = atoi(argv[2]);
	Pr.alpha = atof(argv[3]);
	Pr.beta = atof(argv[4]);
	Pr.gamma = atof(argv[5]);
	Pr.delta = atof(argv[6]);

	char nomefile [80];

	double integration_time;
	vector<int> choice;
	vector<int>::iterator choiceI;
	vector<double> jump_time(DIM_ENSEMBLE);

	vector<double> clfactv(4,1);
	clfactv[0]=exp(-s);
	clfactv[2]=exp(s);

	double next_time;

	vector<Systems> system(DIM_ENSEMBLE);
	choice.reserve(10000000);

	r = gsl_rng_alloc(gsl_rng_taus);

	sprintf(nomefile,"%sASEP%s_%s_%s_%s_%s_%s_%s.txt", argv[8],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7] );
	fout = fopen(nomefile,"w");

	for (i=0; i<DIM_ENSEMBLE; i++) {
	    system[i].part_sum.resize(4);
	    system[i].part_sum.assign(4,0);
/*	    system[i].part_sum_conserved.resize(4);
	    system[i].part_sum_conserved.assign(4,0);
*/	}

	/* * *
	 * START THE SIMULATION
	 */
	weakness = 0;
	steps = 0;

	for (i=0; i<DIM_ENSEMBLE; i++) {
		system[i].current=0;
		system[i].site=0;
/*	    system[i].part_sum[0] = Pr.alpha*(1-system[i].site) *exp(-s);
	    system[i].part_sum[1] = system[i].part_sum[0] + Pr.beta  *system[i].site;
	    system[i].part_sum[2] = system[i].part_sum[1] + Pr.gamma *system[i].site  * exp(s);
	    system[i].part_sum[3] = system[i].part_sum[2] + Pr.delta *(1-system[i].site);
	    system[i].part_sum_conserved[0] = Pr.alpha*(1-system[i].site);
	    system[i].part_sum_conserved[1] = system[i].part_sum_conserved[0] + Pr.beta  * system[i].site;
	    system[i].part_sum_conserved[2] = system[i].part_sum_conserved[1] + Pr.gamma * system[i].site;
	    system[i].part_sum_conserved[3] = system[i].part_sum_conserved[2] + Pr.delta * (1-system[i].site);*/

	    system[i].part_sum[0] = Pr.alpha*(1-system[i].site) ;
	    system[i].part_sum[1] = system[i].part_sum[0] + Pr.beta  *system[i].site;
	    system[i].part_sum[2] = system[i].part_sum[1] + Pr.gamma *system[i].site;
	    system[i].part_sum[3] = system[i].part_sum[2] + Pr.delta *(1-system[i].site);
	}

	for (i=0; i<DIM_ENSEMBLE; i++) {
		jump_time[i] = time_increment(&system[0], i, r);
	}

	cloning = 0;

    i = min_element(jump_time.begin(), jump_time.end()) - jump_time.begin();

	integration_time = T;

	while (jump_time[i] < integration_time){
		steps++;
	    R = gsl_rng_uniform(r) * system[i].part_sum.back();
	    j = upper_bound(system[i].part_sum.begin(),system[i].part_sum.end(),R) - system[i].part_sum.begin();
	    status = update(&system[0], i, j, Pr, s);
	    if ( status != 1 ) {
	    	cerr << "ERROR " << i << " " << j << endl;
	        return EXIT_FAILURE;
	    }
	    next_time = time_increment(&system[0], i, r);
	    jump_time[i] = jump_time[i] + next_time;

	     /* * *
	      * Compute the cloning factor
	      */
	    clfact  = clfactv[j];//system[i].part_sum.back() ;//system[i].part_sum_conserved.back();
	    if (clfact != 1) {
	        y = rand_to_int(clfact,r);
	        if (y> DIM_ENSEMBLE) {
	        	weakness ++;
	        }
	    }
	    else {
	        y=1;
	    }

	    if (y == 0) {
	    	/* * *
	         * Substitution
	         */
	    	// do{num_rand = gsl_rng_uniform_int(r,DIM_ENSEMBLE);}while(num_rand==i);
		    num_rand = gsl_rng_uniform_int(r, DIM_ENSEMBLE - 1);
      		if (num_rand >= i) num_rand++;
	        system[i] = system[num_rand];
	        jump_time[i] = jump_time[num_rand];
	    }
	    else if (y > 1) {
		     /* * *
		      * Cloning and overwriting
		      */
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
	    /* * *
	     * return the label of the next clone that will do a transition.
	     */
	     
	    cloning = cloning + log( (double)(DIM_ENSEMBLE+clfact-1)/(double)DIM_ENSEMBLE);	     
	    i = min_element(jump_time.begin(), jump_time.end()) - jump_time.begin();
	}

    fprintf(fout,"%.2f\t%Le\t%f\t", s, cloning/T, (double)weakness/(double)steps);

	/* * *
	 * Thermodynamic integration
	 */
	mean_observable=0;
	for (int it=0; it<DIM_ENSEMBLE; it++ ) {
		mean_observable = mean_observable + (double)(system[it].current)/T;
	}
	fprintf(fout,"%e\n", mean_observable/(double)DIM_ENSEMBLE);

	gsl_rng_free(r);
	fclose(fout);
	return EXIT_SUCCESS;
}

inline int rand_to_int(double clfact, gsl_rng * r){
    return    (clfact + gsl_rng_uniform(r));
}

