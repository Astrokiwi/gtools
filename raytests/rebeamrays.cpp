#include <omp.h>
#include "prototypes.h"
#include <vector>
#include <list>    
//#include <math.h>
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <deque>
#include <array>

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
    //for ( int ig = 20000 ; ig<22000 ; ig+=1000 ) {
    { int ig = 20000;
        /*if ( ig%100==0 ) {
            std::cout << ig << std::endl;
        }*/
        gasP_t p_cent;
        
        for ( int ii=0 ; ii<3 ; ii++ ) {
            p_cent.r[ii] = 0.;
        }
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
    mintemptree *temptree;
    std::clock_t start;
    double duration;
    double mtot=0.;
    int ignored = 0;

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
    
    // build tree of lists of particles for hitting smoothing lengths
    tree = new octtreelist(rmin,size);


    for ( int ip = 0; ip<d->ng ; ip++ ) {
        mtot+=d->gasP[ip].m;
    }
    
    // only count particles with enough optical depth to have a significant
    // effect, even if *all* of the gas particles had this optical depth
    for ( int ip = 0; ip<d->ng ; ip++ ) {
    //for ( int ip = 0; ip<10 ; ip++ ) {
    //for ( int ip = 0; ip<d->ng ; ip+=1000 ) {
        double max_depth = mtot*kernel->flat_w2(0.,square(d->gasP[ip].h))*dust_physics;
        if ( max_depth>100. ) {
            tree->head_addP(ip,d->gasP[ip].r,d->gasP[ip].h);
        } else {
            ignored++;
        }
    }
    std::cout << "N ignored:" << ignored << std::endl;

    std::cout << "mintemptree:" << std::endl;
    
    // build tree of minimum temperatures of particles for testing
    // if we need to pass rays between them
    double r_cent[3];
    for ( int ii=0 ; ii<3 ; ii++ ) {
        r_cent[ii] = (rmin[ii]+rmax[ii])/2.;
    }
    temptree = new mintemptree(r_cent,size);
    //for ( int ip = 0; ip<50 ; ip++ ) {
    for ( int ip = 0; ip<d->ng ; ip++ ) {
        //d->gasP[ip].uint=sin(ip*100.)*.1+.2;
        temptree->addp(d->gasP[ip].r,ip,d->gasP[ip].uint);
    }
    //temptree->dump();
    
    //std::cout << "DONE!" << std::endl;
    //exit(1);
    
    
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
    //std::array<bool, d->ng> alreadyDone; 
    //{ int ig = 20000;

    start = omp_get_wtime();

    // force into centre and make extra hot
    /*for ( int ii=0 ; ii<3 ; ii++ ) {
        d->gasP[0].r[ii] = 0.;
    }*/
    //d->gasP[0].uint*=1000.;
    //d->gasP[0].uint/=1000;


    #pragma omp parallel default(none) shared(temptree,tree,std::cout,d,lums,fluxes,mtot)
    {
        int nthreads = omp_get_num_threads();
        int ithread = omp_get_thread_num();
        bool alreadyDone[d->ng];
        std::vector< double > local_fluxes ( d->ng,0.);
    //#pragma omp for schedule(dynamic, 1)
        #pragma omp for schedule(dynamic, 1)
        //for ( int ig = 0 ; ig<65536 ; ig++ ) {
        //for ( int ig = 0 ; ig<64 ; ig++ ) {
        for ( int ig = 0 ; ig<d->ng ; ig++ ) {
        //{ int ig = 6400;
        //{ int ig = 0;
            //std::cout << ithread << " " << ig << std::endl;
            if ( ig%100==0 ) {
                std::cout << ig << std::endl;
            }
            gasP_t *p1 = &d->gasP[ig];
        
            // propagate through tree
            // needs physics to make sense
            std::vector<int>* threshList = temptree->getThreshList(d,ig,mtot);
        
            /*if ( ithread==0 ) {
                std::cout << threshList->size() << std::endl;
            }*/
            
                /*threshList = temptree->getThreshList(d,ig,mtot);
                std::cout << threshList->size() << std::endl;
                for ( int ii=0 ; ii<100 ; ii++ ) {
                    std::cout << (*threshList)[ii] << " ";
                }
                std::cout << std::endl;*/
            //}
            
            //exit(1);
        
            //if ( false ) {        
                for ( unsigned int ii=0 ; ii<threshList->size() ; ii++ ) {
                //for ( int jg = 0 ; jg<d->ng ; jg++ ) {
            //for ( int jg = 0 ; jg<d->ng ; jg+=100 ) {
            //for ( int jg = 0 ; jg<d->ng ; jg+=5 ) {
                //for ( int jg = 0 ; jg<10000 ; jg++ ) {
                //{   int jg = 1;
                /*if ( jg%1000==0 ) {
                    std::cout << jg << std::endl;
                }*/
                   int jg = (*threshList)[ii];
                   if ( ig!=jg ) {
                       gasP_t *p2 = &d->gasP[jg];
                       double d12_norm = get_d12_norm(p1,p2);
                       //local_fluxes[jg]+=(*lums)[ig]*exp(-this->one_sph_tree_ray(tree,p1,p2,alreadyDone,d,ig,jg,d12_norm))/d12_norm;
                       //local_fluxes[jg]=this->one_sph_tree_ray(tree,p1,p2,alreadyDone,d,ig,jg,d12_norm);
                       if ( sqrt(d12_norm)<p1->h+p2->h ) {
                        // particles are overlapping - use different method
                        // pretend they aren't touching?
                           local_fluxes[jg]+=   planck_mean_emit_absorb * square(dust_physics) * powerfour(p1->uint)
                                                * p1->m * p2->m
                                                * exp(-this->one_sph_tree_ray(tree,p1,p2,alreadyDone,d,ig,jg,d12_norm))
                                                / (p1->h+p2->h);
                       } else {
                           local_fluxes[jg]+=   planck_mean_emit_absorb * square(dust_physics) * powerfour(p1->uint)
                                                * p1->m * p2->m
                                                * exp(-this->one_sph_tree_ray(tree,p1,p2,alreadyDone,d,ig,jg,d12_norm))
                                                / d12_norm;
                       }
                       //local_fluxes[jg]+=1.;
                   }
                   //fluxes[ig] = 1;
                   //fluxes[jg] = 1;
                   //std::cout << p1->r[0] << " " <<p1->r[1] << " " <<p1->r[2] << " " << std::endl;
                   //std::cout << p2->r[0] << " " <<p2->r[1] << " " <<p2->r[2] << " " << std::endl;
                }
            //}
        
        }   
    
    
     
//     std::cout << "ncomp:" << ncomp << std::endl;
/*        for ( int iproc=0 ; iproc<nthreads ; iproc++ ) {
            #pragma omp for
            for ( int ii=(ithread+iproc)%nthreads ; ii<d->ng ; ii+=nthreads ) {
                fluxes[ii] += local_fluxes[ii];
            }
            #pragma omp barrier
        }*/
        #pragma omp critical
        for ( int ii=0 ; ii<d->ng ; ii++ ) {
            fluxes[ii] += local_fluxes[ii];
        }
            
    }
    duration = ( omp_get_wtime() - start );
    std::cout << "actual tree time:" << duration << std::endl;
    delete tree;
    
    return fluxes;
}

double rebeamrays::one_sph_tree_ray(octtreelist *tree, gasP_t *p1, gasP_t *p2,bool alreadyDone[],gizData_t *d, int ig, int jg, double d12_norm) {
    //double d12_norm = get_d12_norm(p1,p2);
    //SET_TYPE<int >* hitp_set;
    //hitp_set = tree->beam(p1->r,p2->r); // maybe build list of pointers to cells instead of list of particle IDs?
    std::vector<PLIST_TYPE<int>*> *list_of_lists = tree->beam(p1->r,p2->r); // list of pointers to cell particle lists

    /*SET_TYPE<int>::iterator it;
    std::sort(hitp_set->begin(),hitp_set->end());
    it = std::unique(hitp_set->begin(), hitp_set->end());
    hitp_set->resize(std::distance(hitp_set->begin(),it) );*/


    //std::cout << "n to check = " << hitp_set->size() << "/" << d->ng << std::endl;
    //std::vector< bool > alreadyDone ( d->ng,false );
    //std::array< bool > alreadyDone ( d->ng,false );
    //alreadyDone.fill(false);
    //std::unordered_set< int > alreadyDone;

    /*for ( int ig = 0 ; ig<d->ng ; ig++ ) {
        alreadyDone[ig] = false;
    }*/
    memset(alreadyDone, 0, sizeof(bool)*d->ng );
    
    double depth = 0.;
    
    //for ( std::vector<PLIST_TYPE<int>*>::iterator it2 = list_of_lists->begin() ; it2!=list_of_lists->end() ; it2++ ) {
    //    PLIST_TYPE<int > *hitp_set=*it2;
    unsigned int lols = list_of_lists->size();
    //                 std::cout << lols << std::endl;
    for ( unsigned int ilist=0 ; ilist<lols ; ilist++ ) {
        PLIST_TYPE<int > *hitp_set = (*list_of_lists)[ilist];
        // std::cout << ilist << " " << hitp_set << std::endl;
        //for ( SET_TYPE<int >::iterator it = hitp_set->begin() ; it != hitp_set->end(); it++ ) {
        //for ( PLIST_TYPE<int >::iterator it = hitp_set->begin() ; it != hitp_set->end(); it++ ) {
        //    int kg = *it;
            /*if ( kg<0 || kg>=d->ng ) {
                std::cout << "bad kg " << kg << std::endl;
                exit(1);
            }*/
    //                     ncomp++;
        unsigned int hs = hitp_set->size() ;
        for ( unsigned int ii=0 ; ii<hs ; ii++ ) {
            int kg = (*hitp_set)[ii];
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
                        //double invh2 = 1./h2;
       
       
                        if ( d2int<=h2 ) {
                            //fluxes[kg] = 1.;
                            //depth+=1.e-2;
                            //fluxes[kg] = ig;
                            depth+=p3->m*kernel->flat_w2(d2int,h2)*dust_physics;
                            /*if ( depth>100. ) {
                                std::cout<< "BIG DEPTH" << std::endl;
                                std::cout<< ii << " " << p3->m << " " << dust_physics << std::endl;
                                std::cout<< kernel->flat_w2(d2int,h2) << " " << d2int << " " << h2 << std::endl;
                                std::cout<< kernel->flat_w2(0.,h2) << " " << kernel->flat_w2(d2int,1.) << std::endl;
                                exit(1);
                            }*/
                        }
       
                        //fluxes[kg] = d2int;
                    }
                }
            }
            //fluxes[kg] = 3;
        }
    }
    //std::cout << "size: " << alreadyDone.size() << std::endl;
    //fluxes[jg] += (*lums)[ig]/d12_norm*exp(-depth);
    //local_fluxes[jg]+= (*lums)[ig]/d12_norm*exp(-depth);
    //(*local_fluxes)[jg]+= (*lums)[ig]*exp(-depth);
    //fluxes[jg] += (*lums)[ig]*exp(-depth);
    //delete hitp_set;
    delete list_of_lists;
    //exit(1);
    // 
    
    //return (*lums)[ig]*exp(-depth);
    return depth;
}

