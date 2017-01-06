#include <omp.h>
#include "prototypes.h"
#include <vector>
#include <list>    
//#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <deque>

rebeamrays::rebeamrays() {
    kernel = new Kernel();
}

std::vector< double > rebeamrays::SPH_fluxes_thin(gizData_t *d, std::vector< double >* lums) {
    std::vector< double > fluxes ( d->ng,0.);

    for ( int ig = 0 ; ig<d->ng ; ig++ ) {
        for ( int jg = 0 ; jg<d->ng ; jg++ ) {
            if ( ig!=jg ) {
                double rad2 = 0;
                for ( int ii=0 ; ii<3 ; ii++ ) {
                    rad2+=square(d->gasP[ig].r[ii]-d->gasP[jg].r[ii]);
                }
                fluxes[jg]+=(*lums)[ig]/rad2;
            }
        }
    }
    return fluxes;
}


bool rebeamrays::between_particles(gasP_t *p1,gasP_t *p2,gasP_t *p3) {
    //bool between = false;    
    double z1=0.,z2=0.,z3=0.;
    for ( int jj=0 ; jj<3 ; jj++ ) {
        double dd = p1->r[jj]-p2->r[jj];
        z1+=p1->r[jj]*dd;
        z2+=p2->r[jj]*dd;
        z3+=p3->r[jj]*dd;
    }
    if ( z1>z2 ) {
        if ( z1>z3 && z3>z2 ) {
            //between = true;
            return true;
        }
    } else if ( z1<z2 ) {
        if ( z1<z3 && z3<z2 ) {
            //between = true;
            return true;
        }
    } else {
        // ?? should not happen
        std::cout << "z1=z2 ???" << std::endl;
        exit(0);
    }
    //if ( (z1<z2 && z2<z3) || (z1>z2 && z2>z3 ) ) return true;
    //return between;
    return false;
}


double rebeamrays::intersect_d2_nonorm(gasP_t *p1,gasP_t *p2,gasP_t *p3) {
    double dr31[3],dr32[3];
    double crossp;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        dr31[ii] = p3->r[ii]-p1->r[ii];
        dr32[ii] = p3->r[ii]-p2->r[ii];
    }
    
    crossp = square( dr31[1]*dr32[2]-dr31[2]*dr32[1]);
    crossp+= square(-dr31[0]*dr32[2]+dr31[2]*dr32[0]);
    crossp+= square( dr31[0]*dr32[1]-dr31[1]*dr32[0]);
    
    return crossp;
}


double rebeamrays::get_d12_norm(gasP_t *p1,gasP_t *p2) {
    double norm=0.;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        norm+=square(p1->r[ii]-p2->r[ii]);
    }
    
    return norm;
}


std::vector< double > rebeamrays::SPH_fluxes_naive(gizData_t *d, std::vector< double >* lums) {
    std::vector< double > fluxes ( d->ng,0.);

    // check what particles are in the middle
    
    /*int i1 = 700, i2=1000;
    
    gasP_t *p1 = &d->gasP[i1];
    gasP_t *p2 = &d->gasP[i2];
    
    for ( int ig = 0 ; ig<d->ng ; ig++ ) {
        if ( ig!=i1 && ig!=i2 ) {
            gasP_t *p3 = &d->gasP[ig];
            //double dr[3];
            double z1=0.,z2=0.,z3=0.;
            for ( int jj=0 ; jj<3 ; jj++ ) {
                double dd = p1->r[jj]-p2->r[jj];
                z1+=p1->r[jj]*dd;
                z2+=p2->r[jj]*dd;
                z3+=p3->r[jj]*dd;
            }
            if ( z1>z2 ) {
                if ( z1>z3 && z3>z2 ) {
                    fluxes[ig] = 1.;
                }
            } else if ( z1<z2 ) {
                if ( z1<z3 && z3<z2 ) {
                    fluxes[ig] = 1.;
                }
            } else {
                // ?? should not happen
                fluxes[ig] = -1.;
            }
        } else {
            fluxes[ig] = 2.;
        }
    }*/

    // ray from ig to jg passing through kg
    //long ncomp = 0;
    //for ( int ig = 0 ; ig<d->ng ; ig+=10 ) {
    //for ( int ig = 0 ; ig<d->ng ; ig++ ) {
    { int ig = 20000;
    //for ( int ig = 20000 ; ig<22000 ; ig+=1000 ) {
        /*if ( ig%100==0 ) {
            std::cout << ig << std::endl;
        }*/
        gasP_t *p1 = &d->gasP[ig];
        for ( int jg = 0 ; jg<d->ng ; jg++ ) {
        //for ( int jg = 0 ; jg<d->ng ; jg+=100 ) {
        //{   int jg = 1000;
           
           gasP_t *p2 = &d->gasP[jg];
           if ( ig!=jg ) {
               double d12_norm = get_d12_norm(p1,p2);
               double depth = 0;
               for ( int kg = 0 ; kg<d->ng ; kg++ ) {
                   // ncomp++;
                   if ( kg!=ig && kg!=jg ) {
                       // determine if particle is between
                       gasP_t *p3 = &d->gasP[kg];
                       
                       if ( between_particles(p1,p2,p3) ) {
                           double d2int = intersect_d2_nonorm(p1,p2,p3)/d12_norm;
                           
                           double h2 = square(p3->h);
                           
                           if ( d2int<=h2 ) {
                               //fluxes[kg] = ig;
                               depth+=1.e-2;
                           }
                           
                           //fluxes[kg] = d2int;
                       }
                   }
               }
               //fluxes[jg] += (*lums)[ig]/d12_norm*exp(-depth);
               fluxes[jg] += (*lums)[ig]*exp(-depth);
            }
            
           //fluxes[ig] = 1;
           //fluxes[jg] = 1;
        }
    }
//     std::cout << "ncomp:" << ncomp << std::endl;

    return fluxes;
}

std::vector< double > rebeamrays::SPH_fluxes_grid(gizData_t *d, std::vector< double >* lums) {
    std::vector< double > fluxes ( d->ng,0.);
    
    // set up grid
    int L = 32;
    double gridR=.02;
    //double cell_dR = (gridR*2.)/(L);
    //std::vector<std::vector<std::vector<std::list<double> > > > grid;
    gridList *gL = new gridList(L);
    
    std::list<int > *thisList;
    std::list<int >::iterator it;
    std::list<std::list<int >* >::iterator lit;
    
    // place particles in grid
    for ( int ig=0 ; ig<d->ng ; ig++ ) {
        gasP_t *p = &d->gasP[ig];
        int ir[3];
        bool ingrid = true;
        for ( int ii=0 ; ii<3 ; ii++ ) {
            //ir[ii] = floor(((p->r[ii]+gridR)/2./gridR)*L);
            ir[ii] = igridr(p->r[ii],gridR,L);
            if ( ir[ii]<0 || ir[ii]>=L ) {
                ingrid = false;
            }
        }
        if ( ingrid ) {
            // ** TODO - particles get placed in multiple cells based on their smoothing length
            thisList = gL->listAt(ir[0],ir[1],ir[2]);
            thisList->push_back(ig);
        } else {
            // ** TODO - place in "extras" list
        }
    }
    
    /*thisList->push_back(5);
    thisList->push_back(13);
    thisList->push_back(21);
    
    for ( it = thisList->begin() ; it != thisList->end(); it++ ) {
        std::cout << *it << std::endl;
    }*/

    
        
    /*std::ofstream myfile;
    myfile.open("hout.dat");
    for ( int ig=0 ; ig<d->ng ; ig++ ) {
        gasP_t *p1 = &d->gasP[ig];
        myfile << p1->h << " " << ceil(p1->h/cell_dR) << std::endl;
    }
    myfile.close();
    exit(0);    */
    
    // ray from ig to jg passing through kg
    //for ( int ig = 0 ; ig<d->ng ; ig++ ) {
    //for ( int ig = 0 ; ig<d->ng ; ig+=10 ) {
    for ( int ig = 20000 ; ig<22000 ; ig+=1000 ) {
    //{ int ig = 20000;
        /*if ( ig%100==0 ) {
            std::cout << ig << std::endl;
        }*/
        gasP_t *p1 = &d->gasP[ig];
        //for ( int jg = 0 ; jg<d->ng ; jg++ ) {
        //for ( int jg = 0 ; jg<d->ng ; jg+=100 ) {
        {   int jg = 1000;
           gasP_t *p2 = &d->gasP[jg];
           int idest[3];
           for ( int ii=0 ; ii<3 ; ii++ ) {
               idest[ii] = igridr(p2->r[ii],gridR,L);
            }
           if ( ig!=jg ) {
                // move through grid
                std::list< std::list<int >* > cellsHit;
                
                double ddr[3],liner[3],dt[3];
                int iline[3],idir[3];
                int dir_choice;
                bool ingrid=true,hitp=false;
                //int idr[3];
                
                // calculate initial directions
                for ( int ii=0 ; ii<3 ; ii++ ) {
                    ddr[ii] = p2->r[ii]-p1->r[ii];
                    liner[ii] = p1->r[ii];
                    iline[ii] = igridr(p1->r[ii],gridR,L);
                    if ( ddr[ii]>0. ) {
                        idir[ii] = 1;
                    } else if ( ddr[ii]<0. ) {
                        idir[ii] = -1;
                    } else {
                        idir[ii] = 0;
                    }
                }
                
                // add initial cell
                cellsHit.push_back(gL->listAt(iline[0],iline[1],iline[2]));
                
                // loop through line
                while(ingrid && !hitp) {
                    // calculate "dt" to each edge
                    // go in the direction with the lowest "dt"
                    // this will currently fail if ddr[ii]=0. for any direction, so we'll have to fix that
                    for ( int ii=0 ; ii<3 ; ii++ ) {
                        dt[ii] = (rgridi(iline[ii]+idir[ii],gridR,L)-liner[ii])/ddr[ii];
                    }
                    
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
                        cellsHit.push_back(gL->listAt(iline[0],iline[1],iline[2]));
                    }
                }
                
                
                
                
                
                /*
                // just a dumb line
                for ( int ix = 0 ; ix<L ; ix++ ) {
                    cellsHit.push_back(gL->listAt(ix,ix,L/2));
                }*/
                
                
                // find "hit" particles in cells along ray
                double d12_norm = get_d12_norm(p1,p2);
                for ( lit = cellsHit.begin() ; lit != cellsHit.end(); lit++ ) {
                    thisList = (*lit);
                    
                    for ( it = thisList->begin() ; it != thisList->end(); it++ ) {
                        int kg = *it;
                        gasP_t *p3 = &d->gasP[kg];
                        
                        if ( between_particles(p1,p2,p3) ) {
                            double d2int = intersect_d2_nonorm(p1,p2,p3)/d12_norm;
                           
                            double h2 = square(p3->h);
                           
                            if ( d2int<=h2 ) {
                                //fluxes[kg] = 1.;
                                //depth+=1.;
                                fluxes[kg] = 3;
                            }
                           
                            //fluxes[kg] = d2int;
                        }
                        //fluxes[kg] = 3;
                    }
                }
                
                //int ir[3];
           }
           fluxes[ig] = 1;
           fluxes[jg] = 1;
        }
    }
    return fluxes;
}

std::vector< double > rebeamrays::SPH_fluxes_tree(gizData_t *d, std::vector< double >* lums) {
    std::vector< double > fluxes ( d->ng,0.);
    octtreelist *tree;
    std::clock_t start;
    double duration;


    start = omp_get_wtime();
    double rmin[3],rmax[3];
    
    for ( int ii=0 ; ii<3 ; ii++ ) {
        rmin[ii] = d->gasP[0].r[ii];
        rmax[ii] = d->gasP[0].r[ii];
    }

    for ( int ip=1 ; ip<d->ng ; ip++ ) {
        for ( int ii=0 ; ii<3 ; ii++ ) {
            if ( d->gasP[ip].r[ii]<rmin[ii] ) {
                rmin[ii] = d->gasP[ip].r[ii];
            }
            if ( d->gasP[ip].r[ii]>rmax[ii] ) {
                rmax[ii] = d->gasP[ip].r[ii];
            }
        }
    }
    
    double size =rmax[0]-rmin[0];
    for ( int ii=1 ; ii<3 ; ii++ ) {
        size = std::max(size,rmax[ii]-rmin[ii]);
    }



    /*std::cout << "box min/max: "<< std::endl;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        std::cout << rmax[ii] << " " <<rmin[ii] <<std::endl;
    }*/
    
    // box contains [rmin,rmin+size)
    // so must increase size by a small amount so that the "rightmost"
    // particle actually fits strictly inside
    size*=(1.+1.e-6);

    tree = new octtreelist(rmin,size);
    
    for ( int ip = 0; ip<d->ng ; ip++ ) {
    //for ( int ip = 0; ip<10 ; ip++ ) {
    //for ( int ip = 0; ip<d->ng ; ip+=1000 ) {
        tree->head_addP(ip,d->gasP[ip].r,d->gasP[ip].h);
    }
    duration = ( omp_get_wtime() - start );
    std::cout << "treebuildtime:" << duration << std::endl;

    //tree->sort();
    
    //std::cout << d->ng << std::endl;
    
    //tree->matrix_dump();
    //tree->dump();
    //exit(1);
// ray from ig to jg passing through kg
//     long ncomp = 0;
    //for ( int ig = 0 ; ig<d->ng ; ig++ ) {
    //for ( int ig = 0 ; ig<d->ng ; ig+=5000 ) {
    { int ig = 20000;
    //for ( int ig = 20000 ; ig<22000 ; ig+=1000 ) {
        /*if ( ig%100==0 ) {
            std::cout << ig << std::endl;
        }*/
        gasP_t *p1 = &d->gasP[ig];
        //for ( int jg = 0 ; jg<d->ng ; jg+=100 ) {
        for ( int jg = 0 ; jg<d->ng ; jg+=10 ) {
        //{   int jg = 1000;
        if ( jg%1000==0 ) {
            std::cout << jg << std::endl;
        }
           gasP_t *p2 = &d->gasP[jg];
           double d12_norm = get_d12_norm(p1,p2);
           if ( ig!=jg ) {
                SET_TYPE<int >* hitp_set;
                hitp_set = tree->beam(p1->r,p2->r); // maybe build list of pointers to cells instead of list of particle IDs?

                /*SET_TYPE<int>::iterator it;
                std::sort(hitp_set->begin(),hitp_set->end());
                it = std::unique(hitp_set->begin(), hitp_set->end());
                hitp_set->resize(std::distance(hitp_set->begin(),it) );*/
                
                
                //std::cout << "n to check = " << hitp_set->size() << "/" << d->ng << std::endl;
                std::vector< bool > alreadyDone ( d->ng,false );
                //std::unordered_set< int > alreadyDone;
                
                double depth = 0.;

                for ( SET_TYPE<int >::iterator it = hitp_set->begin() ; it != hitp_set->end(); it++ ) {
                    int kg = *it;
                    /*if ( kg<0 || kg>=d->ng ) {
                        std::cout << "bad kg " << kg << std::endl;
                        exit(1);
                    }*/
//                     ncomp++;
                //for ( int ii=0 ; ii<hitp_set->size() ; ii++ ) {
                //    int kg = (*hitp_set)[ii];
                    //if ( !alreadyDone[kg] ) {
                    if ( kg!=ig && kg!=jg ) {
                        if ( !alreadyDone[kg] ) {
                        //if ( alreadyDone.find(kg)==alreadyDone.end() ) {
                            alreadyDone[kg] = true;
                            //alreadyDone.insert(kg);
                            gasP_t *p3 = &d->gasP[kg];
                    
                            if ( between_particles(p1,p2,p3) ) {
                                double d2int = intersect_d2_nonorm(p1,p2,p3)/d12_norm;
                       
                                double h2 = square(p3->h);
                       
                                if ( d2int<=h2 ) {
                                    //fluxes[kg] = 1.;
                                    depth+=1.e-2;
                                    //fluxes[kg] = ig;
                                }
                       
                                //fluxes[kg] = d2int;
                            }
                        }
                    }
                    //fluxes[kg] = 3;
                }
                //std::cout << "size: " << alreadyDone.size() << std::endl;
                fluxes[jg] += (*lums)[ig]/d12_norm*exp(-depth);
                //fluxes[jg] += (*lums)[ig]*exp(-depth);
                delete hitp_set;
                //exit(1);
           }
           //fluxes[ig] = 1;
           //fluxes[jg] = 1;
           //std::cout << p1->r[0] << " " <<p1->r[1] << " " <<p1->r[2] << " " << std::endl;
           //std::cout << p2->r[0] << " " <<p2->r[1] << " " <<p2->r[2] << " " << std::endl;
        }
    }    
//     std::cout << "ncomp:" << ncomp << std::endl;
    
    delete tree;
    
    return fluxes;
}



