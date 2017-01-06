#include "prototypes.h"
#include <vector>
#include <math.h>
#include <iostream>
#include <stdlib.h>

void Kernel::integrate_flat() {
    flattened_table.resize(L_flat);
    std::fill(flattened_table.begin(), flattened_table.end(), 0.);

    flattened_table2.resize(L_flat);
    std::fill(flattened_table2.begin(), flattened_table2.end(), 0.);
    
    double r3;
    double r2;
    double rsq;
    
    double dz = 1./zL;
    
    for ( int ir = 0 ; ir<L_flat ; ir++ ) {
        r2 = ((double)(ir))/L_flat;
        rsq = r2;
        for ( int iz = 0 ; iz<zL ; iz++ ) {
            r3 = sqrt(pow(((double)iz)/zL,2)+pow(r2,2));
            flattened_table[ir]+=2.*w(r3)*dz; // even around 0.

            r3 = sqrt(pow(((double)iz)/zL,2)+rsq);
            flattened_table2[ir]+=2.*w(r3)*dz; // even around 0.
        }
    }
    
    /*for ( int ir = 0 ; ir<L_flat ; ir++ ) {
        std::cout << flattened_table[ir] << std::endl;
    }*/

}

void Kernel::set_L_flat(int L_in, int zL_in) {
    L_flat = L_in;
    zL = zL_in;
    integrate_flat();
}

Kernel::Kernel() {
    //set_L_flat(256,256); // default
    set_L_flat(256,256); // default
}


        
double Kernel::w(double r, double h) {
    double k;
    double x = r/h;
    
    // by default - cubic spline
    if ( x<.5 ) {
        k=1.+6.*(x-1.)*pow(x,2);
    } else if ( x<1. ) {
        k=2.*pow(1.-x,3);
    } else {
        k = 0.;
    }
    k*=8./M_PI/pow(h,3);
    return k;
}


        
double Kernel::w(double x) {
    return w(x,1.);
}

double Kernel::flat_w(double r, double h) {
    // from table
    double x = zL*r/h;
    int ix = (int)(x);
    if ( ix>=zL ) {
        return 0.;
    }
    
    // interpolate
    double klow,khigh,kf;
    klow = flattened_table[ix];
    kf = x-ix;
    if ( ix==zL-1 ) {
        khigh = 0.;
    } else {
        khigh = flattened_table[ix+1];
    }
    
    double k = klow*(1.-kf)+khigh*kf;
    /*if ( k<0. ) {
        std::cout << klow << " " << khigh << " " << kf << std::endl;
        exit(1);
    }*/
    
    k/=square(h); // flattened kernel has h=1, need to take real h into account
    return k;
}



double Kernel::flat_w2(double r2, double h2) {
    // from table
    double x = zL*r2/h2;
    int ix = (int)(x);
    if ( ix>=zL ) {
        return 0.;
    }
    
    // interpolate
    double klow,khigh,kf;
    klow = flattened_table2[ix];
    kf = x-ix;
    if ( ix==zL-1 ) {
        khigh = 0.;
    } else {
        khigh = flattened_table2[ix+1];
    }
    
    double k = klow*(1.-kf)+khigh*kf;
    /*if ( k<0. ) {
        std::cout << klow << " " << khigh << " " << kf << std::endl;
        exit(1);
    }*/
    
    k/=h2; // flattened kernel has h=1, need to take real h into account
    return k;
}



double Kernel::flat_w2_q(double r2, double invh2) {
    // from table
    double x = zL*r2*invh2;
    int ix = (int)(x);
    if ( ix>=zL ) {
        return 0.;
    }
    
    // interpolate
    double klow,khigh,kf;
    klow = flattened_table2[ix];
    kf = x-ix;
    if ( ix==zL-1 ) {
        khigh = 0.;
    } else {
        khigh = flattened_table2[ix+1];
    }
    
    //double k = klow+(khigh-klow)*kf;
    /*if ( k<0. ) {
        std::cout << klow << " " << khigh << " " << kf << std::endl;
        exit(1);
    }*/
    
    //k*=invh2; // flattened kernel has h=1, need to take real h into account
    return klow+(khigh-klow)*kf*invh2;
}
