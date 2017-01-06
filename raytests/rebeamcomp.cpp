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
            myfile << (*depths)[ig] << std::endl;
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
    }

    
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
    fluxes = dr->SPH_fluxes_tree(myData,&lums);
    duration = ( omp_get_wtime() - start );
    std::cout<<"tree version:  "<< duration <<'\n';
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
