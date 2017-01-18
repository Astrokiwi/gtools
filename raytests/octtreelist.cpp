#include <omp.h>
#include "prototypes.h"
#include <list>    
#include <unordered_set>
#include <set>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string> 
#include <limits>
#include <deque>

/* Oct-tree for finding collisions between a ray and SPH particles.
Particles added to the tree are propagated down the oct-tree nodes
(creating new nodes along the way if possible) until they reach the 
shallowest level of refinement where cells are larger than the diameter
of the particle, i.e. 2*h. The particle is added to the particle-lists
for all cells at that refinement level that it intersects with, which
is up to eight cells.

This means that particles can be placed in branch nodes as well as in
leaf nodes. This stops particles with large softening lengths from being
placed in a very large number of cells - the max is 8.

A beam can be propagated through the tree. The beam is propagated through
all cells at each level of refinement in tern, returning a list of pointers
to particle-lists of cells intersected by the beam. This list will contain
duplicates.
*/

// create head node
octtreelist::octtreelist(double r_left[3], double size) {
    for ( int ii=0 ; ii<3 ; ii++ ) {
        this->r_left[ii] = r_left[ii];
        this->ixyz[ii] = 0;
    }
    this->size = size;
    this->level = 0;
    //myP = new PLIST_TYPE<int>;
    //myP.reserve(100);
    //std::cout << "head size, coords=" << size << " " << r_left[0] << " " << r_left[1] << " " << r_left[2] << std::endl;
    //myP.reserve(1e4);

    this->nodes = new octtreelist*[8];
    for ( int ii=0 ; ii<8 ; ii++ ) {
        this->nodes[ii] = NULL;
    }
    
    this->maxlevel = 0;

}

// create leaf node
octtreelist::octtreelist(int level, int ixyz[3]) {
    for ( int ii=0 ; ii<3 ; ii++ ) {
        this->ixyz[ii] = ixyz[ii];
    }
    this->level = level;
    //myP = new PLIST_TYPE<int>;
    //myP.reserve(100);
    //myP.reserve(1e4);
    
    this->nodes = new octtreelist*[8];
    for ( int ii=0 ; ii<8 ; ii++ ) {
        this->nodes[ii] = NULL;
    }

}

void octtreelist::sort() {
    //myP->sort();
    for ( int ii=0 ; ii<8 ; ii++ ) {
        if ( nodes[ii] ) {
            nodes[ii]->sort();
        }
    }
}

void octtreelist::dump() {
    for ( int ii=0; ii<this->level; ii++ ) {
        std::cout << ".";
    }
    //std::cout << myP->size() << std::endl;
    std::cout << myP.size() << std::endl;
    for ( int ii=0; ii<8; ii++ ) {
        if ( nodes[ii] ) {
            nodes[ii]->dump();
        } else {
            for ( int ii=0; ii<this->level; ii++ ) {
                std::cout << ".";
            }
            std::cout << "-" << std::endl;
        }
    }
}

void octtreelist::matrix_dump() {
    for ( int ilevel = 0 ; ilevel<=this->maxlevel ; ilevel++ ) {
        std::ofstream f;
        f.open("treedump"+std::to_string(ilevel)+".dat");
        int L = 1<<ilevel;
        for ( int ix=0 ; ix<L ; ix++ ) {
            for ( int iy = 0 ; iy<L ; iy++ ) {
                int sizeTot = 0;
                for ( int iz=0 ; iz<L ; iz++ ) { 
                    octtreelist* thisCell = this->getCellAt(ilevel,ix,iy,iz);
                    if ( !thisCell ) {
                        sizeTot+=thisCell->getPList()->size();
                    }
                }
                f << sizeTot << " ";
            }
            f << std::endl;
        }
        f.close();
    }
}

PLIST_TYPE<int>* octtreelist::getPList() {
    //return myP;
    // std::cout << "myP:" << &myP << std::endl;
//     std::cout << "this:" << this << std::endl;
    return &(this->myP);
}

octtreelist* octtreelist::getCellAt(int ilevel,int ix, int iy, int iz) {
    int ixyz[3];
    ixyz[0] = ix;
    ixyz[1] = iy;
    ixyz[2] = iz;
    return this->getCellAt(ilevel,ixyz);
}
    
octtreelist* octtreelist::getCellAt(int ilevel,int ixyz[3]) {
    if ( ilevel==this->level ) {
//         std::cout << "cellat: " << ilevel << " " << this << std::endl;
        return this;
    } else {
        // determine which node
        int inode,dlevel;
        int ixyz_rel[3];
        /*ixyz[0] = ix;
        ixyz[1] = iy;
        ixyz[2] = iz;*/
        dlevel = ilevel-this->level-1;
        for ( int ii=0 ; ii<3 ; ii++ ) {
            ixyz_rel[ii] = (ixyz[ii]>>dlevel)-(this->ixyz[ii]<<1);
        }
        //std::cout << std::endl;
        inode=ixyz_rel[0]+(ixyz_rel[1]<<1)+(ixyz_rel[2]<<2);

        if ( inode>=8 || inode<0 ) {
            return NULL; // out of bounds
        }
        // check if node exists
        if ( !nodes[inode] ) {
            return NULL; // no cell there
        }
        return nodes[inode]->getCellAt(ilevel,ixyz);
    }
}

void octtreelist::addPCoord(int ip,int ilevel,int ix,int iy,int iz) {
    /*if ( ip==1954 ) {
        std::cout << ip << " " << ilevel << " " << ix << " " << iy << " " << iz << std::endl;
    }*/
    //if ( this->level==0 ) {
        //std::cout << ip << " " << this->level << " " << ilevel << " " << ix << " " << iy << " " << iz << std::endl;
        //std::cout << ip << " " << ilevel << " " << this->level << " " << this->ixyz[0] << " " << this->ixyz[1] << " " << this->ixyz[2] << std::endl;
    //}
    if ( ilevel==this->level ) {
        // add to myself on this level
        //myP->push_back(ip);
        myP.push_back(ip);
        //std::cout << "nodevector: " << myP.capacity() << " " << myP.size() << std::endl;
    } else {
        // add to subnode
        
        // determine which node
        int inode,dlevel;
        int ixyz[3],ixyz_rel[3];
        ixyz[0] = ix;
        ixyz[1] = iy;
        ixyz[2] = iz;
        dlevel = ilevel-this->level-1;
        /*std::cout << " incoords:";
        for ( int ii=0 ; ii<3 ; ii++ ) {
            std::cout << ixyz[ii] << " ";
        }
        std::cout << std::endl;
        std::cout << " mycoords:";
        for ( int ii=0 ; ii<3 ; ii++ ) {
            std::cout << this->ixyz[ii] << " ";
        }
        std::cout << std::endl;
        std::cout << "relcoords:";*/
        for ( int ii=0 ; ii<3 ; ii++ ) {
            //ixyz_rel[ii] = ixyz[ii]-(this->ixyz[ii]>>dlevel);
            ixyz_rel[ii] = (ixyz[ii]>>dlevel)-(this->ixyz[ii]<<1);
            //std::cout << ixyz_rel[ii] << " ";
        }
        //std::cout << std::endl;
        // ixyz_rel should be =0 or =1 for all three entries
        inode=ixyz_rel[0]+ixyz_rel[1]*2+ixyz_rel[2]*4;
        //std::cout << "-->add to node:" << inode << std::endl;
        
        if ( inode>=8 || inode<0 ) {
            if ( this->level==0 ) {
                // This just means the cell is out of bounds.
                // This happens if the smoothing length overlaps the
                // edge of the tree of cells. If the tree contains all 
                // particles, then we can ignore this, because no ray
                // will pass outside this box when going between two
                // particles
                return;
            }
            // if ilevel!=0 then something is wrong
            std::cout << ip << " -->add to node:" << inode << std::endl;
            std::cout << ilevel << " " << this->level << " " << dlevel << std::endl;
            for ( int ii=0 ; ii<3 ; ii++ ) {
                std::cout << ixyz_rel[ii] << " " << ixyz[ii] << " " << this->ixyz[ii] << " " << (ixyz[ii]>>dlevel) << " " << (this->ixyz[ii]<<1) << std::endl;
            }
            exit(1);
            return;
        }

        /*std::cout << "inode: " << inode << std::endl;
        for ( int jnode=0 ; jnode<8 ; jnode++ ) {
            //std::cout << (nodes[jnode]==NULL) << " ";
            std::cout << nodes[jnode] << " ";
        }
        std::cout << std::endl;
        std::cout << &myP << std::endl;*/

        // check if node exists
        if ( !nodes[inode] ) {
            //std::cout << "addcoords:" << std::endl;
            int ixyz_node[3];
            for ( int ii=0 ; ii<3 ; ii++ ) {
                ixyz_node[ii] = (this->ixyz[ii]<<1)+ixyz_rel[ii];
                //std::cout << ixyz_node[ii] << " " << (this->ixyz[ii]<<1) << " " <<ixyz_rel[ii] <<std::endl;
            }
            nodes[inode] = new octtreelist(this->level+1,ixyz_node);
            //std::cout << "NEW NODE " << inode << std::endl;
        }
        //std::cout << nodes[inode] << std::endl;
        nodes[inode]->addPCoord(ip,ilevel,ix,iy,iz);
    }
}

// should only be called on the head node
// *requires* that particle is in grid
void octtreelist::head_addP(int ip, double r_p[3], double h) {
    int depth_cap = 8; // max is 256^2
    //int depth_cap = 4;

    // goes in level where size of cell<hs
    int p_level;
    p_level = (int)floor(log2(size/h))-1; // +1 because we want 2h within each cell
    if ( p_level>depth_cap ) {
        p_level = depth_cap;
    }
    //std::cout << "p_level=" << p_level << " h:" << h << std::endl;
    /*if ( ip==1954 ) {
        std::cout << ip << " " << r_p[0] << " " << r_p[1] << " " << r_p[2] << std::endl;
        std::cout << r_left[0]+size << " " << r_left[1]+size << " " << r_left[2]+size << std::endl;
    }*/
    
    this->maxlevel = std::max(p_level,this->maxlevel);
    
    if ( p_level<=0 ) {
        // add to myself, the top level
        //myP->push_back(ip);
        myP.push_back(ip);
        //std::cout << "headvector: " << myP.capacity() << " " << myP.size() << std::endl;
    } else {
        // add to lower levels
        // find which cells and level this particle should go
        double lsize = this->size/(1<<p_level); // cell size at target level
        int cent_coords[3]; // integer coordinates of left corner of cell at target level containing particle centre
        double diff_coords[3]; // displacement from corner of centre cell to particle centre
        int bump_over[3]; // which direction (if any) the particle overlaps with adjacent cels
        //std::cout << "lsize=" << lsize << " " << (1<<p_level) << std::endl;
        //std::cout << (1>>2) << " " << (1<<2) << std::endl;
        //std::cout << "coords, cent_coords =";

        for ( int ii=0 ; ii<3 ; ii++ ) {
            cent_coords[ii] = (int)floor((r_p[ii]-r_left[ii])/lsize);
            //std::cout << r_p[ii] << " " << cent_coords[ii] << std::endl;
        }
        // add to centre cell
        this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1],cent_coords[2]);
        //exit(1);
        // add to adjacent cells if particle overlaps. Loop over all 8 adjacent cells
        for ( int ii=0 ; ii<3 ; ii++ ) {
            diff_coords[ii] = r_p[ii]-r_left[ii]-lsize*cent_coords[ii];
            //diff_coords[ii] = r_p[ii]-rtreei(cent_coords[ii],ii,p_level);
            if ( h+diff_coords[ii]>=lsize ) {
                bump_over[ii] = 1;
            } else if ( diff_coords[ii]-h<0. ) {
                bump_over[ii] = -1;
            } else {
                bump_over[ii] = 0;
            }
            //std::cout << diff_coords[ii] << " " << diff_coords[ii]/lsize << " " << h/lsize << " " << bump_over[ii] << std::endl;
        }
        //exit(1);
        //return;
        // place in all adjacent cells
        if ( bump_over[0]!=0 ) {
            this->addPCoord(ip,p_level,cent_coords[0]+bump_over[0],cent_coords[1],cent_coords[2]);
            if ( bump_over[1]!=0 ) {
                this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1]+bump_over[1],cent_coords[2]);
                this->addPCoord(ip,p_level,cent_coords[0]+bump_over[0],cent_coords[1]+bump_over[1],cent_coords[2]);
                if ( bump_over[2]!=0 ) {
                    this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1]+bump_over[1],cent_coords[2]+bump_over[2]);
                    this->addPCoord(ip,p_level,cent_coords[0]+bump_over[0],cent_coords[1],cent_coords[2]+bump_over[2]);
                    this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1],cent_coords[2]+bump_over[2]);
                    this->addPCoord(ip,p_level,cent_coords[0]+bump_over[0],cent_coords[1]+bump_over[1],cent_coords[2]+bump_over[2]);
                }
            }
        } else if ( bump_over[1]!=0 ) {
            this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1]+bump_over[1],cent_coords[2]);
            if ( bump_over[2]!=0 ) {
                this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1]+bump_over[1],cent_coords[2]+bump_over[2]);
                this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1],cent_coords[2]+bump_over[2]);
            }
        } else if ( bump_over[2]!=0 ) {
            this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1],cent_coords[2]+bump_over[2]);
        }
    }
}

// return the set of all particles in cells that intersect the beam
// only called on head
std::vector<PLIST_TYPE<int>*>* octtreelist::beam(double r1[3], double r2[3]) {
//SET_TYPE<int > octtreelist::beam(double r1[3], double r2[3]) {
    //SET_TYPE<int > *phitset = new SET_TYPE<int >();
    //SET_TYPE<int > *phitset_in;
    std::vector<PLIST_TYPE<int>*> *list_of_lists = new std::vector<PLIST_TYPE<int>*>;
    for ( int il=0 ; il<=maxlevel ; il++ ) {
        beam_level(il,r1,r2,list_of_lists);

        //phitset_in = beam_level(il,r1,r2);
        // merge unique values
        //phitset->merge(phitset_in);
        /*for ( std::unordered_set<int >::iterator it = phitset_in->begin() ; it!=phitset_in->end(); it++ ) {
            if ( phitset->find(*it)==phitset_in->end() ) {
                phitset
            }
        }*/
        //phitset->insert(phitset_in->begin(),phitset_in->end());
        //phitset->merge(*phitset_in);
        //phitset->insert(phitset->end(), phitset_in->begin(), phitset_in->end());
        //list_of_lists->insert(phitset->end(), phitset_in->begin(), phitset_in->end());
        //delete phitset_in;
    }
    //return phitset;
    return list_of_lists;
}


// find grid coordinate at this level
int octtreelist::itreer(double r, int idim, int ilevel) {
    return (int)floor(((r-this->r_left[idim])/this->size)*(1<<ilevel));
}

// find "left" corner of cell at given coordinate
double octtreelist::rtreei(int ii, int idim, int ilevel) {
    return (ii*this->size)/(1<<ilevel)+this->r_left[idim];
}

// return the set of all particles in cells at this particular level of refinement that intersect the beam
// only called on head
//SET_TYPE<int >* octtreelist::beam_level(int ilevel, double r1[3], double r2[3]) {
//    SET_TYPE<int > *phitset_out = new SET_TYPE<int >();
void octtreelist::beam_level(int ilevel, double r1[3], double r2[3],std::vector<PLIST_TYPE<int>*> *list_of_lists) {
    
   //int ncells =0;

   int idest[3];
   for ( int ii=0 ; ii<3 ; ii++ ) {
       idest[ii] = itreer(r2[ii],ii,ilevel);
    }

    
    // Move through cells at this level
    double ddr[3],liner[3],dt[3];
    int iline[3],idir[3];
    int dir_choice;
    bool ingrid=true,hitp=false;
    
    int L = 1<<ilevel;
    
    // calculate initial directions
    for ( int ii=0 ; ii<3 ; ii++ ) {
        ddr[ii] = r2[ii]-r1[ii];
        liner[ii] = r1[ii];
        iline[ii] = itreer(r1[ii],ii,ilevel);
        if ( ddr[ii]>0. ) {
            idir[ii] = 1;
        } else if ( ddr[ii]<0. ) {
            idir[ii] = -1;
        } else {
            idir[ii] = 0;
        }
    }
    
    // add initial cell
    octtreelist *cell0 = this->getCellAt(ilevel,iline[0],iline[1],iline[2]);
    if ( cell0 ) {
        //phitset_out->insert(cell0->getPList()->begin(),cell0->getPList()->end());
        //phitset_out->merge(*(cell0->getPList()));
        //phitset_out->insert(phitset_out->end(), cell0->getPList()->begin(),cell0->getPList()->end());
        // std::cout << ilevel << " " << cell0 << std::endl;
        list_of_lists->push_back(cell0->getPList());
    }
    //ncells = 1;
    
    // Loop through line
    // This version is inefficient and loops through lots of empty cells, especially at greatest refinement level!
    // TODO: make less dumb
    while(ingrid && !hitp) {
        //std::cout << ilevel << " " << iline[0] << " " <<iline[1] << " " << iline[2] << " " << phitset_out->size() << std::endl;
        //std::cout << " " << liner[0] << " " << liner[1] << " " << liner[2] << " " << std::endl;
        //std::cout << " " << iline[0] << " " << iline[1] << " " << iline[2] << " " << std::endl;
        //std::cout << " " << idir[0] << " " << idir[1] << " " << idir[2] << " " << std::endl;
        // calculate "dt" to each edge
        // go in the direction with the lowest "dt"
        // this will currently fail if ddr[ii]=0. for any direction, so we'll have to fix that
        for ( int ii=0 ; ii<3 ; ii++ ) {
            if ( idir[ii]==1 )         dt[ii] = (rtreei(iline[ii]+1,ii,ilevel)-liner[ii])/ddr[ii]; // "right" side of this cell is "left" cell of next cell
            else if ( idir[ii]==-1 )   dt[ii] = (rtreei(iline[ii],ii,ilevel)-liner[ii])/ddr[ii]; // "left" side of this cell
            else dt[ii] = std::numeric_limits<double>::max();//(rtreei(iline[ii]+idir[ii],ii,ilevel)-liner[ii])/ddr[ii];
            //std::cout << rtreei(iline[ii]+idir[ii],ii,ilevel);
        }
        //std::cout << " " << idir[0] << " " << idir[1] << " " << idir[2] << " " << std::endl;
        //std::cout << " " << dt[0] << " " << dt[1] << " " << dt[2] << " " << std::endl;
        //std::cout << std::endl;
        dir_choice = 0;
        for ( int ii=1 ; ii<3 ; ii++ ) {
            if ( dt[ii]<dt[dir_choice] ) {
                dir_choice = ii;
            }
        }
        
        // jump to the next cell
        iline[dir_choice]+=idir[dir_choice];
        
        // update line location
        for ( int ii=0 ; ii<3 ; ii++ ) {
            liner[ii]+=dt[dir_choice]*ddr[ii];
        }
        
        // check for finishing conditions
        ingrid = true;
        hitp = true;
        for ( int ii=0 ; ii<3 ; ii++ ) {
            if ( iline[ii]<0 || iline[ii]>=L ) ingrid=false;
            if ( iline[ii]!=idest[ii] ) hitp = false;
        }
        
        // add the list of particles
        if ( ingrid ) {
            //cellsHit.push_back(gL->listAt(iline[0],iline[1],iline[2]));
            octtreelist *hitCell = this->getCellAt(ilevel,iline[0],iline[1],iline[2]);
            if ( hitCell ) {
                //phitset_out->insert(hitCell->getPList()->begin(),hitCell->getPList()->end());
                //phitset_out->merge(*(hitCell->getPList()));
                //phitset_out->insert(phitset_out->end(), hitCell->getPList()->begin(),hitCell->getPList()->end());
                //std::cout << ilevel << " " << hitCell << std::endl;
                list_of_lists->push_back(hitCell->getPList());
            }
            //ncells++;
        }
    }
    //std::cout << "level:" << ilevel << " ncells:" << ncells << std::endl;
    //return phitset_out;
}







