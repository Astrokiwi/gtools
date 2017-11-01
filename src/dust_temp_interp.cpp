#include <stdlib.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "dust_temp_interp.h"

/*// helper functions
double square(double x) {
    return x*x;
}


int square(int x) {
    return x*x;
}*/

/*! find index in (short) table where value fits to the right of
return -1 if it's to the left of everything, and n-1 if it's to the right of everything
*/
int CoolHeatTab::value_to_index(double edges[], int nedges, double value) {
    int i;
    
    if ( value<edges[0] ) return -1;
    for ( i=1 ; i<nedges ; i++ ) {
        if ( value<edges[i] ) {
            return i-1;
        }
    }
    return nedges-1;
}


/*! weighting for interpolation
No error checking!
*/
double CoolHeatTab::right_weight(double edges[], int nedges, int index, double value) {
    if ( index<0 ) return 1.;
    if ( index>=nedges-1) return 0.;
    
    return (value-edges[index])/(edges[index+1]-edges[index]);
}

int CoolHeatTab::agn_tab_index(int id, int it, int ii, int is) {
    return  is +
            agn_nsurf*ii +
            agn_nsurf*agn_nintensity*it +
            agn_nsurf*agn_nintensity*agn_ntemp*id;
};
// 
// struct[] coolHeatDust CoolHeatTab::interpTabVec(double density[], double temperature[], double intensity[], double coldens[], int n) {
//     struct outp[n];
//     
//     for ( int ii=0 ; ii<n ; ii++ ) {
//         outp[ii] = interpTab(density[ii],temperature[ii],intensity[ii],coldens[ii]);
//     }
//     
//     return outp;
// } 

// takes input values in *cgs units*
struct coolHeatDust CoolHeatTab::interpTab(double density, double temperature, double intensity, double coldens) {
    int jj,kk;
    
    int tab_index[4];
    int interp_indices[4][2];
    //int interp_right_index[4];
    double interp_weights[4][2];
    double values[4];
    
    struct coolHeatDust outp;
    
    values[0] = log10(density);
    values[1] = log10(temperature);
    values[2] = log10(intensity);
    values[3] = coldens;
    
    for (jj=0 ; jj<4 ; jj++ ) {
        tab_index[jj] = value_to_index(tables[jj],*ntabs[jj],values[jj]);
        interp_indices[jj][0] = std::max(tab_index[jj],0);
        interp_indices[jj][1] = std::min(tab_index[jj]+1,*ntabs[jj]-1);
        interp_weights[jj][1] = right_weight(tables[jj],*ntabs[jj],tab_index[jj],values[jj]);
        interp_weights[jj][0] = 1.-interp_weights[jj][1];
    }

    // too clever by half
    double heat_interp=0.,cool_interp=0.,dust_interp=0.,dg_interp=0.;
    double weightsum = 0.;
    for (jj=0 ; jj<16 ; jj++ ) {
        int idex = agn_tab_index(interp_indices[0][jj%2],interp_indices[1][(jj/2)%2],interp_indices[2][(jj/4)%2],interp_indices[3][(jj/8)%2]);
        //int idex = agn_tab_index(5,5,interp_indices[2][(jj/4)%2],5);
        double weight = 1.;
        for ( kk=0 ; kk<4 ; kk++ ) {
            weight*=interp_weights[kk][(jj/(1<<kk))%2];
        }
        //weight = interp_weights[2][(jj/4)%2];
        
        heat_interp+=weight*agn_heat_tab[idex];
        cool_interp+=weight*agn_cool_tab[idex];
        dust_interp+=weight*agn_dust_tab[idex];
        dg_interp+=weight*agn_dg_tab[idex];
        weightsum+=weight;
//         std::cout << weight << " " << idex << " " << agn_dust_tab[idex] << " " << agn_cool_tab[idex] << " " << agn_heat_tab[idex] << " " << dust_interp << std::endl;

    }
    
    outp.dCool = cool_interp;
    outp.dHeat = heat_interp;
    outp.dustT = dust_interp;
    outp.dg = dg_interp;
    
    outp.id = interp_indices[0][0];
    outp.fd = interp_weights[0][1];
    outp.it = interp_indices[1][0];
    outp.ft = interp_weights[1][1];
    outp.ii = interp_indices[2][0];
    outp.fi = interp_weights[2][1];
    outp.ic = interp_indices[3][0];
    outp.fc = interp_weights[3][1];
    
    
//     if ( outp.dustT>1.e9 ) {
//         std::cout << "bad dustT" << std::endl;
//         exit(1);
//     }

//     if ( outp.dg>1. ) {
//         std::cout << "bad dg" << std::endl;
//         exit(1);
//     }
    
    return outp;
}

//CoolHeatTab::CoolHeatTab(std::string flabels,std::string ftab) {
CoolHeatTab::CoolHeatTab(const char* flabels,const char* ftab) {
    int i,ntab;
    
    double arad_filler;
    
    // read in tables
    //std::ifstream f_tablab (flabels.c_str(), std::ifstream::in);
    std::ifstream f_tablab (flabels, std::ifstream::in);

    f_tablab >> agn_ndense;
    agn_dense_vals = new double[agn_ndense];
    for ( i=0 ; i<agn_ndense ; i++ ) {
        f_tablab >> agn_dense_vals[i];
    }

    f_tablab >> agn_ntemp;
    agn_temp_vals = new double[agn_ntemp];
    for ( i=0 ; i<agn_ntemp ; i++ ) {
        f_tablab >> agn_temp_vals[i];
    }

    f_tablab >> agn_nintensity;
    agn_intensity_vals = new double[agn_nintensity];
    for ( i=0 ; i<agn_nintensity ; i++ ) {
        f_tablab >> agn_intensity_vals[i];
    }

    f_tablab >> agn_nsurf;
    agn_surf_vals = new double[agn_nsurf];
    for ( i=0 ; i<agn_nsurf ; i++ ) {
        f_tablab >> agn_surf_vals[i];
    }
    
    f_tablab.close();
    
    ntab=agn_ndense*agn_ntemp*agn_nintensity*agn_nsurf;
    agn_heat_tab = new double[ntab];
    agn_cool_tab = new double[ntab];
    agn_dust_tab = new double[ntab];
    agn_dg_tab = new double[ntab];
    
    //std::ifstream f_tabvals (ftab.c_str(), std::ifstream::in);
    std::ifstream f_tabvals (ftab, std::ifstream::in);
    
    for ( i=0 ; i<ntab ; i++ ) {
        f_tabvals >> agn_heat_tab[i] >> agn_cool_tab[i] >> agn_dust_tab[i] >> arad_filler >> agn_dg_tab[i];
//         if ( agn_heat_tab[i]==0 ) {
//             std::cout << "borked" << i << std::endl;
//             exit(0);
//         }
//         if ( agn_dg_tab[i]>1. ) {
//             std::cout << "borked" << i << std::endl;
//             exit(0);
//         }
    }
    f_tabvals.close();
    
    //** Unit conversions - keep in cgs!
    
    // 
//     // number density in cm**-3 to physical density in internal units
//     const double molecular_mass = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */
//     double logdenseconvert = log10(
//                                     (molecular_mass*PROTONMASS)      // from cm**-3 to g/cm**3
//                                     /All.UnitDensity_in_cgs
//                                     );
//     for ( i=0 ; i<agn_ndense ; i++ ) {
//         agn_dense_vals[i]+=logdenseconvert;
//     }
//     
//     // temperature in K to internal energy in internal units
//     double logtempconvert = log10(All.UnitMass_in_g / All.UnitEnergy_in_cgs*BOLTZMANN/GAMMA_MINUS1/(molecular_mass*PROTONMASS));
//     for ( i=0 ; i<agn_ntemp ; i++ ) {
//         agn_temp_vals[i]+=logtempconvert;
//     }
// 
//     // intensity in erg/s/cm**2 to internal units
//     double logintensityconvert = log10(square(All.UnitLength_in_cm)/All.UnitEnergy_in_cgs*All.UnitTime_in_Megayears*SEC_PER_MEGAYEAR);
//     for ( i=0 ; i<agn_nintensity ; i++ ) {
//         agn_intensity_vals[i]+=logintensityconvert;
//     }
//     
//     // surface density in cm**-2 to mass surface density in internal units
//     double surfconvert = (molecular_mass*PROTONMASS)    // from cm**-2 to g/cm**-2
//                           /All.UnitMass_in_g
//                           *square(All.UnitLength_in_cm);
//     for ( i=0 ; i<agn_nsurf ; i++ ) {
//         agn_surf_vals[i]*=surfconvert;
//     }
//     
//     // heating/cooling from erg/s/cm**3 to internal units energy/time/volume
//     double heatcoolconvert = log10(pow(All.UnitLength_in_cm,3)/All.UnitEnergy_in_cgs*All.UnitTime_in_Megayears*SEC_PER_MEGAYEAR);
//     for ( i=0 ; i<ntab ; i++ ) {
//         agn_heat_tab[i]+=heatcoolconvert;
//         agn_cool_tab[i]+=heatcoolconvert;
//     }
//     
    
    // set up pointers for better looping
    tables[0] = agn_dense_vals;
    ntabs[0] = &agn_ndense;
    tables[1] = agn_temp_vals;
    ntabs[1] = &agn_ntemp;
    tables[2] = agn_intensity_vals;
    ntabs[2] = &agn_nintensity;
    tables[3] = agn_surf_vals;
    ntabs[3] = &agn_nsurf;

}



// int main(int argc, char **argv) {
//     CoolHeatTab *oldTable = new CoolHeatTab("coolheat_tab_marta/shrunk_table_labels_050517.dat","coolheat_tab_marta/shrunk_table_050517.dat");
//     struct coolHeatDust outp;
//     
//     outp = oldTable->interpTab(1.e-1,1.e7,1.e8,1.e16);
//     
//     std::cout << "Heat:" << outp.dHeat << " Cool:" << outp.dCool << " Dust:" << outp.dustT << std::endl;
// }

