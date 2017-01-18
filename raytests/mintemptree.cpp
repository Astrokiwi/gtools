#include "prototypes.h"

#include <algorithm>
#include <iostream>

/* Oct-tree of minimum temperatures for determining which particles
receive non-negligible heating. Here, particles are placed in unique
leaf nodes.
*/


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
        double halfsize = size/2;
        for ( int ii=0 ; ii<3 ; ii++ ) {
            if ( r_p[ii]>this->r_cent[ii] ) {
                r_cent[ii] = this->r_cent[ii]+halfsize;
            } else {
                r_cent[ii] = this->r_cent[ii]-halfsize;
            }
        }
        this->nodes[inode] = new mintemptree(r_cent,halfsize);
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
