#include <omp.h>
#include "prototypes.h"
#include <vector>
#include <list>    
#include <math.h>
#include <stdlib.h>

gridList::gridList(int L) {
    this->L = L;
    listVector = new std::vector<std::list< int>>(L*L*L);
}

std::list<int >* gridList::listAt(int ix, int iy, int iz) {
    return &((*listVector)[ix+iy*L+iz*L*L]);
}
