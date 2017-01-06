#include "hdf5.h"
#include <stdio.h>
#include <math.h>
#include <fstream>

#include "prototypes.h"


int main() {
    gizmo_hdf5_reader *gread = new gizmo_hdf5_reader("/export/1/djw/gizmo_public/disc_nocool_midres_out/snapshot_096.hdf5");

    gread->readGas();
    
    gizData_t *myData = gread->getData();
    
    /*for ( int ig=0 ; ig<myData->ng ; ig+=5000 ) {
        printf("%g %g %g %g %g %g\n",
                            log10(myData->gasP[ig].rho),
                            log10(myData->gasP[ig].uint),
                            log10(myData->gasP[ig].m),
                            myData->gasP[ig].r[0],
                            myData->gasP[ig].r[1],
                            myData->gasP[ig].r[2]
                            );
    }*/
    
    // Surface density plot
    int L = 512;
    double sd_width = 2.e-2;
    double sd[L][L];
    
    for ( int ix=0 ; ix<L ; ix++ ) {
        for ( int iy = 0 ; iy<L ; iy++ ) {
            sd[ix][iy] = 0.;
        }
    }
    
    for ( int ig=0 ; ig<myData->ng ; ig++ ) {
        gasP_t *myP = &myData->gasP[ig];
        
        int ixy[2];
        bool onscreen = true;
        for ( int jj=0 ; jj<2 ; jj++ ) {
            ixy[jj] = (myP->r[jj]/sd_width+.5)*L;
            if ( ixy[jj]<0 || ixy[jj]>=L ) {
                onscreen = false;
            }
        }
        
        if ( onscreen ) {
            //sd[ixy[0]][ixy[1]]+=myP->m;
            sd[ixy[0]][ixy[1]]=myP->rho;
        }
    }

    std::ofstream myfile;
    myfile.open("matrix.txt");
    for ( int ix=0 ; ix<L ; ix++ ) {
        for ( int iy = 0 ; iy<L ; iy++ ) {
            myfile << sd[ix][iy] << " ";
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open("dots.txt");
    for ( int ig=0 ; ig<myData->ng ; ig++ ) {
        gasP_t *myP = &myData->gasP[ig];
        myfile << myP->r[0] << " " << myP->r[1] << " " << myP->rho<< " " << myP->m << std::endl;
    }
    myfile.close();


    return 1;
}
