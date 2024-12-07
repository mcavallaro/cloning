#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "type_defs.h"
#include "indexed_heap.hpp"



int update(Systems*, int , int, int);
int assign_part_sum(Systems *,int );
inline int rand_to_int(double, gsl_rng*);


/*int reservoir_sampling(gsl_rng * r,std::vector<int> &choices, int N, int k ) {
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
}*/





//this function implements the WRS, Algorithm A-Res of [Efraimidis, Spirakis,  Information Processing Letters 97, (2006) 181 ]
//with Heap
int heap_weighted_reservoir_sampling(gsl_rng * r, std::vector<int> &choices, const std::vector<double> &weights, int k ) {
    int N = weights.size();
    if (k > N)
        return EXIT_FAILURE;
    Heapselect keys(k);
    double key;
    int i;

    for (i = 0; i < k; i++){
        keys.add_large( pow(gsl_rng_uniform(r) , 1./weights[i]), i ) ;
    }    
    // Iterate from the (k+1)th element to nth element
    for (; i < N; i++) {
        key = pow(gsl_rng_uniform(r) ,1./weights[i]);
        if (key > keys.smaller_value(0) ) {
            keys.add_large(key,i);
        }
    }
    std::copy(keys.indexes.begin(), keys.indexes.end(), choices.begin());
    return EXIT_SUCCESS;
}


void usage(char *argv[]){
  printf("********************************************************************\n"
         "**                                                                **\n"
         "**     -*-   Ring Discrete Time Cloning -*-                       **\n"
         "**                                                                **\n"
         "********************************************************************\n"
         "\n\n"
         " This is Free Software - You can use and distribute it under \n"
         " the terms of the GNU General Public License, version 3 or later\n\n"
         " (c) Massimo Cavallaro (m.cavallaro@qmul.ac.uk)\n"
         " Warning: now alpha, beta, gamma, delta are probabilities---alpha+delta=beta+gamma=1\n");
  printf("Usage: %s [# simulation time] [DIM. ENSEMBLE] [N] [s] [output file prx] [output folder]\n\n" , argv[0]);
}

int main(int argc, char *argv[]) {
    if (argc != 7){
        usage(argv);
        fprintf(stderr,"%d\n",argc);
        exit(EXIT_FAILURE);
    }

//     double  mean_observable;


    FILE* filep;

    char file_name[120];

    int j;
    int t=0;


    int i; //span the ensemble
    double R;
    int DIM_ENSEMBLE;
    double T;

    int N = atoi(argv[3]);

    double s = atof(argv[4]);
    double cloning;
    double  clfact;
    int status = 0;
    gsl_rng * r;


    //int y,  num_rand, new_copies;
    char nome [10];

    snprintf(nome,10,"%s",argv[5]);


    vector<int> choice;
    vector<int>::iterator choiceI;

    snprintf(file_name, 120, "%s%s_%s_%s_%s_%s.txt",  argv[6], nome,argv[1], argv[2], argv[3], argv[4]);
    filep = fopen(file_name, "w");

    T = atof(argv[1]);
    DIM_ENSEMBLE = atoi(argv[2]);


    choice.reserve(DIM_ENSEMBLE);
    vector<Systems> sistema(DIM_ENSEMBLE);
    vector<Systems> sistema_tmp(DIM_ENSEMBLE);


    r = gsl_rng_alloc(gsl_rng_taus);

    for (i =0; i< DIM_ENSEMBLE; i++) {
        sistema[i].position=1;
        sistema[i].current=0;
        sistema[i].part_sum.reserve(2);
        sistema[i].part_sum.resize(2);
        sistema[i].part_sum.assign(2,0);
    }


    vector<double> clfvect(DIM_ENSEMBLE);

    cloning = 0;
    t=0;
    while(t<T){
        t++;
        for(i=0; i<DIM_ENSEMBLE; i++){
            status = assign_part_sum(&sistema[0],i);
            if (status != EXIT_SUCCESS){
                fprintf(stderr,"PROBLEM in assign_part_sum %d", i);
                return EXIT_FAILURE;
            }

            R = gsl_rng_uniform(r) * sistema[i].part_sum.back();
            j = upper_bound(sistema[i].part_sum.begin(),sistema[i].part_sum.end(),R) - sistema[i].part_sum.begin();
            status = update(&sistema[0], i, j, N);
            if ( status != 1 ) {
                cerr << "ERROR " << i << " " << j << endl;
                return EXIT_FAILURE;
            }


            if (j){
	            clfvect[i] = exp(s);
            }
            else {
            	clfvect[i] = exp(-s);
            }
//            cout << clfvect[i] << " " ;
        }
//        y = rand_to_int(clfact,r);
        clfact = accumulate(clfvect.begin(),clfvect.end(),0.0);
//        cout <<  clfact << endl;



        status = heap_weighted_reservoir_sampling(r, choice, clfvect, DIM_ENSEMBLE );
        if (status != EXIT_SUCCESS){
            std::cerr << "WRS failed" << std::endl;
            return EXIT_FAILURE;
        }

        for(i=0; i<DIM_ENSEMBLE; i++){
            sistema_tmp[i] = sistema[choice[i]];
        }
        sistema = sistema_tmp;

/*        if (clfact != 1) {
            y = rand_to_int(clfact,r);
        }
        else {
            y=1;
        }
        if (y == 0) {  
            do{num_rand = gsl_rng_uniform_int(r,DIM_ENSEMBLE);}while(num_rand==i);
            sistema[i] = sistema[num_rand];
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
                    sistema[*choiceI] = sistema[i];
                }
            }
        }*/

//        cloning = cloning + log( (double)(DIM_ENSEMBLE+clfact-1)/(double)DIM_ENSEMBLE);
        cloning = cloning + log( (double)(clfact)/(double)DIM_ENSEMBLE);
    }

    fprintf(filep,"%.2f\t%e\n", s, cloning/T);

   /* Thermodynamic integration*/
/*    mean_observable=0;
    for (int it=0; it<DIM_ENSEMBLE; it++ ) {
        mean_observable = mean_observable + (double)(sistema[it].current)/T;
    }
    fprintf(filep,"%e\n", mean_observable/(double)DIM_ENSEMBLE);*/


    fclose(filep);
    gsl_rng_free(r);
    return EXIT_SUCCESS;


}




inline int rand_to_int(double clfact, gsl_rng * r){
    return    (clfact + gsl_rng_uniform(r));
}

