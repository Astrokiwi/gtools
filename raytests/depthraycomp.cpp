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
    
        gasP_t *p1 = &data->gasP[ig];
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
    int tid,nthreads;
/*#pragma omp parallel private(tid,nthreads)
    {

    tid = omp_get_thread_num();
    printf("Hello World from thread = %d\n", tid);
    if (tid == 0) 
    {
        nthreads = omp_get_num_threads();
        printf("Number of threads = %d\n", nthreads);
    }
    }
    exit(1);*/

    //gizmo_hdf5_reader *gread = new gizmo_hdf5_reader("/export/1/djw/gizmo_public/disc_nocool_midres_out/snapshot_096.hdf5");
    gizmo_hdf5_reader *gread = new gizmo_hdf5_reader("/export/1/djw/gizmo_public/disc_nocool_out/snapshot_010.hdf5");

    gread->readGas();
    
    gizData_t *myData = gread->getData();
    
    // make 2D for viz
    /*for ( int ig=0 ; ig<myData->ng ; ig++ ) {
        myData->gasP[ig].r[2] = 0.;
        // add wall for lulz
        if ( ig%2==0 ) {
            myData->gasP[ig].r[0]=.004;
            myData->gasP[ig].r[1]/=2.;
        }
    }*/
    
    depthrays *dr = new depthrays();
    
    std::clock_t start;
    double duration;
    std::vector< double > depths;


    start = omp_get_wtime();
    depths = dr->SPH_depth2(myData);
    duration = ( omp_get_wtime() - start );
    std::cout<<"sqr  version: "<< duration <<'\n';
    savePoints("depth1.dat",&depths,myData);
    
    //exit(1);

    start = omp_get_wtime();
    depths = dr->SPH_depth(myData);
    duration = ( omp_get_wtime() - start );
    std::cout<<"sqrt version: "<< duration <<'\n';
    savePoints("depth2.dat",&depths,myData);
}
