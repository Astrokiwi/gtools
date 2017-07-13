#include "prototypes.h"

#include <algorithm>
#include <iostream>
#include <vector>

/* Oct-tree of minimum temperatures for determining which particles
receive non-negligible heating. Here, particles are placed in unique
leaf nodes.
*/


// a good value for this depends on the particular simulation and is not at all a universal value!
// units of inverse mass (1/(1e10 Msun))
//static const double lum_cutoff = 1.e7;
//static const double lum_cutoff = 5.e4;
//static const double lum_cutoff = 1.e3;
static const double lum_cutoff = 1000.;


// create head node
mintemptree::mintemptree(double r_cent[3], double size) {
    for ( int ii=0 ; ii<3 ; ii++ ) {
        this->r_cent[ii] = r_cent[ii];
    }
    this->size = size;
    
    this->nodes = new mintemptree*[8];
    for ( int ii=0 ; ii<8 ; ii++ ) {
        this->nodes[ii] = NULL;
    }
    
    this->ip = -1; // empty leaf node
    this->minTemp = -1;
}

void mintemptree::propagateThreshList(gizData_t* d, int ip,std::vector<int>* threshList, int ilevel, double mtot) {
    // Does this node fit the criterion?
    // If so - add my particle to the list if I'm a leaf node, or open up my child nodes if I'm a branch node
    
    // calculate distance
    double halfsize = this->size/2;
    double r2 = 0.;
    
    if ( this->ip<=-1 ) {
        // distance to closest edge of cell
        for ( int ii=0 ; ii<3 ; ii++ ) {
            r2+=square(std::max(std::abs(d->gasP[ip].r[ii]-this->r_cent[ii])-halfsize,0.));
            //std::cout << d->gasP[ip].r[ii] << "-" << this->r_cent[ii] << " ";
        }
        //std::cout << std::endl;
    } else {
        // distance to centre of particle I contain
        for ( int ii=0 ; ii<3 ; ii++ ) {
            r2+=square(d->gasP[ip].r[ii]-d->gasP[this->ip].r[ii]);
            //std::cout << d->gasP[ip].r[ii] << "-" << this->r_cent[ii] << " ";
        }
    }
    
    double my_crit;
    
    if ( r2>0. ) {
        my_crit=mtot*dust_physics*planck_mean_emit_absorb/planck_mean_absorb*powerfour(d->gasP[ip].uint/this->minTemp)/r2;
    } else {
        my_crit = lum_cutoff+1.;
    }
    
//     for ( int ii=0 ; ii<ilevel ; ii++ ) {
//         std::cout << ".";
//     }
//     std::cout << my_crit << "/" << lum_cutoff << " " << this->minTemp << " " << r2 << " " << this->size << " " << this->ip << std::endl;
    
    if ( my_crit > lum_cutoff ) {
        if ( this->ip<=-1 ) {
            // propagate to lower nodes
            for ( int ii=0 ; ii<8 ; ii++ ) {
                if ( this->nodes[ii] != NULL ) {
                    this->nodes[ii]->propagateThreshList(d,ip,threshList,ilevel+1,mtot);
                }
            }
        } else {
            threshList->push_back(this->ip);
            /*for ( unsigned int ii=0 ; ii<threshList->size() ; ii++ ) {
                std::cout << (*threshList)[ii] << " ";
            }*/
            //std::cout << std::endl;
        }
    }
    return;
}

std::vector<int>* mintemptree::getThreshList(gizData_t* d, int ip, double mtot) {
    std::vector<int>* threshList = new std::vector<int>;
    
    this->propagateThreshList(d,ip,threshList,0,mtot);
    
    return threshList;
}

void mintemptree::dump() {
    this->dump(0);
}

void mintemptree::dump(int ilevel) {
    for ( int ii=0 ; ii<ilevel ; ii++ ) {
        std::cout << ".";
    }
    std::cout << this->minTemp << " " << this->ip << std::endl;
    for ( int inode=0 ; inode<8 ; inode++ ) {
        if ( this->nodes[inode]!=NULL ) {
            this->nodes[inode]->dump(ilevel+1);
        }
    }
}


void mintemptree::placeInChildNode(double r_p[3], int ip, double temp) {
    int inode=0;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        if ( r_p[ii]>this->r_cent[ii] ) {
            inode+=1<<ii;
        }
    }
    
    if ( this->nodes[inode]==NULL ) {
        double r_cent[3];
        double quartersize = size/4;
        for ( int ii=0 ; ii<3 ; ii++ ) {
            if ( r_p[ii]>this->r_cent[ii] ) {
                r_cent[ii] = this->r_cent[ii]+quartersize;
            } else {
                r_cent[ii] = this->r_cent[ii]-quartersize;
            }
        }
        this->nodes[inode] = new mintemptree(r_cent,size/2.);
    }
    
    this->nodes[inode]->addp(r_p,ip,temp);
}

void mintemptree::addp(double r_p[3], int ip, double temp) {
    if ( this->ip>=0 ) {
        // convert from leaf node to branch node
        
        this->placeInChildNode(this->r_p,this->ip, this->minTemp);
        this->ip = -2; // branch node
        
        this->placeInChildNode(r_p,ip,temp);
        this->minTemp = std::min(temp,this->minTemp);
    } else if ( this->ip==-1 ) {
        // add a particle to empty leaf node
        this->ip = ip;
        for ( int ii=0 ; ii<3 ; ii++ ) {
            this->r_p[ii] = r_p[ii];
        }
        this->minTemp = temp;
    } else {
        // branch node
        this->placeInChildNode(r_p,ip,temp);
        this->minTemp = std::min(temp,this->minTemp);
    }
}
