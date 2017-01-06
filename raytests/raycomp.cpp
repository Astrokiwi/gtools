#include "hdf5.h"
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string>

#include "prototypes.h"

void saveMatrix(std::string fname, std::vector< std::vector<double> >* g) {
    std::ofstream myfile;
    myfile.open(fname.c_str());
    int L = g->size(); // assume square
    for ( int ix=0 ; ix<L ; ix++ ) {
        for ( int iy = 0 ; iy<L ; iy++ ) {
            myfile << (*g)[ix][iy] << " ";
        }
        myfile << std::endl;
    }
    myfile.close();
}

int main() {
    gizmo_hdf5_reader *gread = new gizmo_hdf5_reader("/export/1/djw/gizmo_public/disc_nocool_midres_out/snapshot_096.hdf5");

    gread->readGas();
    
    gizData_t *myData = gread->getData();
    
    raytracers *rt = new raytracers();
    
    rt->setL(1024);
    rt->setw(2.e-2);
    
    std::vector< std::vector<double> > g;
    
    g = rt->naive_surf(myData);
    saveMatrix("sumtest.dat",&g);

    std::vector< std::vector<double> > gsph1;
    gsph1 = rt->SPH_surf(myData);
    saveMatrix("SPHtest.dat",&gsph1);

    std::vector< std::vector<double> > gsph2;
    gsph2 = rt->SPH_surf2(myData);
    saveMatrix("SPHtest2.dat",&gsph2);

    std::vector< std::vector<double> > gsph3;
    gsph3 = rt->SPH_surf3(myData);
    saveMatrix("SPHtest3.dat",&gsph3);

    
    int L = gsph2.size();
    std::vector< std::vector<double> > dg ( L, std::vector<double> (L,0) );
    for ( int ix=0 ; ix<L ; ix++ ) {
        for ( int iy=0 ; iy<L ; iy++ ) {
            dg[ix][iy]=(gsph3[ix][iy]-gsph2[ix][iy])/gsph3[ix][iy];
        }
    }
    saveMatrix("SPHdelta.dat",&dg);

}
