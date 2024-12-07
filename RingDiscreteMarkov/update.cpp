#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cstdio>
#include <iomanip>
#include "type_defs.h"


int assign_part_sum(Systems *system,int i){

    system[i].part_sum[0] = 0.6;
    system[i].part_sum[1] = system[i].part_sum[0] + 0.4;

    return EXIT_SUCCESS;
}


int update(Systems *system, int i, int j, int N){
    int next;
    switch (j){
        case 0:
            system[i].position = (system[i].position+1)%N;
            system[i].current++;            
            break;
        case 1:
            next = system[i].position-1;
            if (next < 0){
                next = N-1;
            }
            system[i].position = next;
            system[i].current--;
            break;
        default:
            fprintf(stderr,"uffa %d\n",j);
            exit(EXIT_FAILURE);
    }

	return 1;
}

