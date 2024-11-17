#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "type_defs.h"

int update(Systems *system, int i, int j, Probabilities Pr, double s){
	if (j == 0) { // alpha
		system[i].current++;		
		system[i].site++ ;
	}
	else if (j == 1){ // beta
		system[i].site--;
	}
	else if (j == 2){ // gamma
		system[i].current-- ;
		system[i].site--;
	}
	else if (j == 3){ // delta
		system[i].site++;
	}
	else {
        fprintf(stderr,"error 777\n" );
        exit(EXIT_FAILURE);
    }

/*    system[i].part_sum[0] = Pr.alpha * (1-system[i].site) * exp(-s);
    system[i].part_sum[1] = system[i].part_sum[0] + Pr.beta  * system[i].site ;
    system[i].part_sum[2] = system[i].part_sum[1] + Pr.gamma * system[i].site * exp(s);
    system[i].part_sum[3] = system[i].part_sum[2] + Pr.delta *(1-system[i].site);
    system[i].part_sum_conserved[0] = Pr.alpha*(1-system[i].site);
    system[i].part_sum_conserved[1] = system[i].part_sum_conserved[0] + Pr.beta  * system[i].site;
    system[i].part_sum_conserved[2] = system[i].part_sum_conserved[1] + Pr.gamma * system[i].site;
    system[i].part_sum_conserved[3] = system[i].part_sum_conserved[2] + Pr.delta * (1-system[i].site);*/


    system[i].part_sum[0] = Pr.alpha * (1-system[i].site) ;
    system[i].part_sum[1] = system[i].part_sum[0] + Pr.beta  * system[i].site ;
    system[i].part_sum[2] = system[i].part_sum[1] + Pr.gamma * system[i].site ;
    system[i].part_sum[3] = system[i].part_sum[2] + Pr.delta *(1-system[i].site);


	return 1;
}