#include <cstdlib>
#include <vector>
#include <cstdio>
#include "type_defs.h"
#include "functions.h"


int update(Systems *system, int i, int j, Parameters br, int dim_chain){
    if (j == 0){
    	system[i].node[0]++;
    	system[i].N++;
    	system[i].current++;
//        system[i].history = system[i].total_time;  
    }
    else if (j < dim_chain){
    	system[i].node[j-1]--;
    	system[i].node[j]++;
    }
    else if (j == dim_chain){
    	system[i].N--;
        system[i].node[j-1]--;
    }
    else {
		fprintf(stderr,"error_777\n");
		exit(EXIT_FAILURE);
	}
	return EXIT_SUCCESS;
}


int assign_part_sum(Systems *system, int i, Parameters br, int dim_chain, vector<double> &times){
    int j;
    if(system[i].node[0]==1){
        system[i].part_sum[0] = 0;
    }
    else if (system[i].node[0]==0){
        system[i].part_sum[0] = w_0(&system[0],i, br, times);
    }
    else{
        fprintf(stderr,"error n_particles is %d\n", system[i].node[0]);
        exit(EXIT_FAILURE);
    }

//    printf( "dd %d %f %d\n",0, system[i].part_sum[0], system[i].node[0]);

    for ( j=1; j<dim_chain; j++) {
        system[i].part_sum[j] = system[i].part_sum[j-1] + system[i].node[j-1] *(1-system[i].node[j]);

//        printf( "%d %f %d\n",j, system[i].part_sum[j],  system[i].node[j]);

    }
    system[i].part_sum[j] = system[i].part_sum[j-1] +  system[i].node[j-1]*br.beta;
//    printf( "%d %f\n",j, system[i].part_sum[j]);

    return EXIT_SUCCESS;
}










