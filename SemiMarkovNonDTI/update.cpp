#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "type_defs.h"
#include "time_jumps.h"

int assign_part_sum(Systems *system,int i){

    system[i].part_sum[0] = Palpha(&system[0],i) * system[i].site;
    system[i].part_sum[1] = system[i].part_sum[0] + Pbeta(&system[0],i)  * (1 - system[i].site);
    system[i].part_sum[2] = system[i].part_sum[1] + Pgamma(&system[0],i) *    system[i].site;
    system[i].part_sum[3] = system[i].part_sum[2] + Pdelta(&system[0],i) * (1 - system[i].site);

    return EXIT_SUCCESS;
}


int update(Systems *system, int i, int j){
	if (j == 0) { // alpha
		system[i].site--;
	}
	else if (j == 1){ // beta
		system[i].site++;
	}
	else if (j == 2){
		system[i].current++;
		system[i].site--;
	}
	else if (j == 3){
		system[i].current--;
		system[i].site++;
	}
	else {
        fprintf(stderr,"error 777\n" );
        exit(EXIT_FAILURE);
    }
	return 1;
}


