#include <omp.h>
#include "hdf5.h"
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <ctime>
#include <stdlib.h>

#include "prototypes.h"

void savePoints(std::string fname, std::vector< double >* depths, gizData_t* data) {
    std::ofstream myfile;
    myfile.open(fname.c_str());
    for ( int ig=0 ; ig<data->ng ; ig++ ) {
    
        //gasP_t *p1 = &data->gasP[ig];
        /*if (    p1->r[0]>0. && p1->r[0]<0.01 &&
                p1->r[1]>0. && p1->r[1]<0.01 ) {*/
            for ( int jj = 0 ; jj<3 ; jj++ ) {
                myfile << data->gasP[ig].r[jj] << " ";
            }
            myfile << (*depths)[ig] << " " << data->gasP[ig].h << std::endl;
        //}
    }
    myfile.close();
}

int main() {
    //int tid,nthreads;
    gizmo_hdf5_reader *gread = new gizmo_hdf5_reader("/export/1/djw/gizmo_public/disc_nocool_midres_out/snapshot_096.hdf5");
    //gizmo_hdf5_reader *gread = new gizmo_hdf5_reader("/export/1/djw/gizmo_public/disc_nocool_out/snapshot_010.hdf5");

    gread->readGas();
    
    gizData_t *myData = gread->getData();
    
    // make 2D for viz
    for ( int ig=0 ; ig<myData->ng ; ig++ ) {
        //myData->gasP[ig].r[2]*=1.e-6;
        //myData->gasP[ig].r[2] = sin(ig)*.01;
        //myData->gasP[ig].r[2] = 0.;
        //myData->gasP[ig].h = 1.e-3;
        // add wall for lulz
        /*if ( ig%3==0 ) {
            myData->gasP[ig].r[0]=.004;
            myData->gasP[ig].r[1]/=2.;
        }*/
        //myData->gasP[ig].r[0]+=1.e-3*sin(myData->gasP[ig].r[0]*3.e2);
        //myData->gasP[ig].r[1]+=1.e-3*sin(myData->gasP[ig].r[1]*3.e2);
    }
    
    // make some artificial clumps
    /*for ( int ig=0 ; ig<myData->ng ; ig+=(myData->ng/20) ) {
        gasP_t *p1 = &myData->gasP[ig];
        for ( int ii=0 ; ii<3 ; ii++ ) {
            p1->r[ii]/=1.5;
        }
        for ( int jg=0 ; jg<myData->ng ; jg++ ) {
            gasP_t *p2 = &myData->gasP[jg];
            double r2=0;
            for ( int ii=0 ; ii<3 ; ii++ ) {
                r2+=square(p2->r[ii]-p1->r[ii]);
            }
            double weight = std::min(1.e-6/r2,.5);
            for ( int ii=0 ; ii<3 ; ii++ ) {
                p2->r[ii]+=(p1->r[ii]-p2->r[ii])*weight;
            }
            
        }
    }*/
    
    //std::cout << "Mass in 10^10 Msun:" << myData->gasP[0].m << std::endl;
    //exit(1);

    
    rebeamrays *dr = new rebeamrays();
    
    std::clock_t start;
    double duration;
    std::vector< double > fluxes; // received fluxes
    std::vector< double > lums(myData->ng); // emitted luminosities
    int nhot = 0;
    for ( int ig=0 ; ig<myData->ng ; ig++ ) {
        lums[ig] = 1.e6;
        /*if ( (ig)%(myData->ng/10)==0 ) {
            lums[ig] = 1.;
            nhot++;
        } else {
            lums[ig] = 1.e-3;
        }*/
    }
    std::cout<<"nhot:"<<nhot<<std::endl;

    start = omp_get_wtime();
    fluxes = dr->SPH_fluxes_tree(myData,&lums); // fluxes are actually dE/dt
    duration = ( omp_get_wtime() - start );
    std::cout<<"tree version:  "<< duration <<'\n';
    // apply cooling
    
    for ( int ig=0 ; ig<myData->ng ; ig++ ) {
        gasP_t *myP = &(myData->gasP[ig]);
        fluxes[ig]+= - myP->m * dust_physics * planck_mean_absorb * powerfour(myP->uint);
    }
    
    savePoints("fluxes3.dat",&fluxes,myData);

    /*start = omp_get_wtime();
    fluxes = dr->SPH_fluxes_grid(myData,&lums);
    duration = ( omp_get_wtime() - start );
    std::cout<<"grid version:  "<< duration <<'\n';
    savePoints("fluxes2.dat",&fluxes,myData);*/
    
    /*start = omp_get_wtime();
    fluxes = dr->SPH_fluxes_naive(myData,&lums);
    duration = ( omp_get_wtime() - start );
    std::cout<<"naive version: "<< duration <<'\n';
    savePoints("fluxes1.dat",&fluxes,myData);*/
    
}
