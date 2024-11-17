
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include "type_defs.h"
#include "functions.h"
#include "update.h"
#include "time_increments_f.h"

int reservoir_sampling(gsl_rng * r, vector<int> &choices, int N, int k ) {
    if (k > N) return EXIT_FAILURE;
    int j,i;
    choices.clear();
    for (i = 0; i < k; i++){
        choices.push_back(i);
    }
    for (; i < N; i++)
    {
        j = gsl_rng_uniform_int(r, i+1);
        if (j < k) {
          choices[j] = i;
        }
    }
    return EXIT_SUCCESS;
}

inline int rand_to_int(double clfact, gsl_rng * r){
    return    (clfact + gsl_rng_uniform(r));
}

void usage(char *argv[]){
  printf("********************************************************************\n"
         "**                                                                **\n"
         "**         -*-            open  ASEP CLONING         -*-          **\n"
         "**                                                                **\n"
         "********************************************************************\n"
         "\n\n"
         " This is Free Software - You can use and distribute it under \n"
         " the terms of the GNU General Public License, version 3 or later\n\n"
         " (c) Massimo Cavallaro (m.cavallaro@qmul.ac.uk)\n\n");
  printf("Ustimes: %s [# simulation time] [DIM. ENSEMBLE] [alpha0] [beta] [dim_chain] [s] [output folder]\n\n" , argv[0]);
}

int main(int argc, char *argv[]) {

  if (argc != 8){
      usage(argv);
      exit(EXIT_FAILURE);
  }

  int dim_chain ;

  double mean_current;
  double * clfact;
  double cloning = 0;
  double s;
  double R;
  double num_clones =0;
  double steps = 0;

  int new_copies;
  FILE* file;

  Parameters br;

  char file_name[120];

  int i;
  int j;
  int y, num_rand;

  int T, DIM_ENSEMBLE;
  int status = 0;
  gsl_rng * r;


  vector<int> choice;
  vector<int>::iterator choiceI;

  choice.reserve(10000000);

  snprintf(file_name, 120, "%smem_TASEP_clon_%s_%s_%s_%s_%s_%s.txt", argv[7],argv[1],argv[2], argv[3], argv[4], argv[5], argv[6]);
  file = fopen(file_name, "w");

  T = atof(argv[1]);
  DIM_ENSEMBLE = atoi(argv[2]);

  vector<double> densityprofile(DIM_ENSEMBLE,0);


  br.alpha0 = atof(argv[3]);
  br.a     = 0.1;
  br.beta  = atof(argv[4]);

  dim_chain = atof(argv[5]);

  s = atof(argv[6]);

  clfact = new double [dim_chain+1];

  for (i =1; i< dim_chain+1; i++) {
    clfact[i] = 1;
  }
  clfact[0] = exp(-s);

  vector<Systems> system(DIM_ENSEMBLE);
  vector<double> times(DIM_ENSEMBLE);

  r = gsl_rng_alloc(gsl_rng_taus);

  for (i =0; i< DIM_ENSEMBLE; i++) {
      system[i].node.reserve(dim_chain);
      system[i].node.resize(dim_chain);
      system[i].age = 0;
      system[i].current=0;
      system[i].part_sum.reserve(dim_chain+1);
      system[i].part_sum.resize(dim_chain+1);
      system[i].part_sum.assign(dim_chain+1,0);
  }

  for (i =0; i< DIM_ENSEMBLE; i++) {
    for (j=0;j<dim_chain;j++){
      if (gsl_rng_uniform(r)>0.5){
        system[i].node[j]=1;
        system[i].N++;
      }
      else 
        system[i].node[j]=0;
    }
  }

  for  (i =0; i< DIM_ENSEMBLE; i++) {
      times[i] = 1;
  }

  for  (i =0; i< DIM_ENSEMBLE; i++) {
      system[i].age  = time_increment_f(&system[0], i, br, r,  times);
      times[i] = 1 + system[i].age ; // time_increment_f(&system[0], i, br, r,  times);
  }



  i = min_element(times.begin(), times.end()) - times.begin();
  while (times[i]  < T ) {
//    cout << times[i] << endl;
    status = assign_part_sum(&system[0],i,br, dim_chain, times);
    if (status != EXIT_SUCCESS){
      fprintf(stderr,"PROBLEM in assign_part_sum %d", i);
      return EXIT_FAILURE;
    }
    R = gsl_rng_uniform(r) * system[i].part_sum.back();
    j = upper_bound(system[i].part_sum.begin(),system[i].part_sum.end(), R)-system[i].part_sum.begin();

    status = update(&system[0], i,j, br, dim_chain);
    if (status!=EXIT_SUCCESS) {
      fprintf(stderr,"ERROR %d %f %d", i, times[i], j);
      return EXIT_FAILURE;
    }

    system[i].age =  time_increment_f(&system[0], i, br, r,  times); //depends on history = $t + t_0$
    times[i] = times[i] + system[i].age;
    
    if (clfact[j] != 1) {
        y = rand_to_int(clfact[j],r);
    }
    else {
        y=1;
    }


    if (y == 0) {
      // do{
      //   num_rand = gsl_rng_uniform_int(r,DIM_ENSEMBLE);
      // }while(num_rand==i);
      num_rand = gsl_rng_uniform_int(r, DIM_ENSEMBLE - 1)
      if (num_rand >= i) num_rand++
      system[i] = system[num_rand];
      times[i] = times[num_rand];
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
          times[*choiceI]=times[i];
        }
      }
    }

    num_clones =num_clones + y;
    steps++;

    cloning = cloning + log( (DIM_ENSEMBLE+clfact[j] -1)/(double)DIM_ENSEMBLE);
    i = min_element(times.begin(), times.end()) - times.begin();
  }
  fprintf(file,"%.2f\t%e\t", s, cloning/T );

 /* * *
  * Thermodynamic integration
  */
  mean_current = 0;
  for (i =0; i< DIM_ENSEMBLE; i++){
    mean_current = mean_current + system[i].current;
  }
  fprintf(file, "%f\t%f\n", mean_current/(double)DIM_ENSEMBLE/(double)T , num_clones/steps);

  delete [] clfact;

  fclose(file);
  gsl_rng_free(r);

}
