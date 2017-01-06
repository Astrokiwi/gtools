#include <omp.h>
#include "prototypes.h"
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <fstream>

depthrays::depthrays() {
    kernel = new Kernel();
}

// CONTINUE DEBUGGERING FAST VERSION
// it's not the big bits that are making the difference
// it's the tail end bits of lots of little particles
// how am I missing these? how bad is it?

std::vector< double > depthrays::SPH_depth(gizData_t *d) {
    std::vector< double > dep ( d->ng,0.);
    
    // "scatter" form, do naive N^2 sum
    for ( int ig = 0 ; ig<d->ng ; ig++ ) {
    //for ( int ig = 0 ; ig<1 ; ig++ ) {
        if ( ig%1000==0 ) std::cout << ig << std::endl;
//         std::ofstream myfile;
//         myfile.open("contrib1.dat");
//    for ( int ig = 0 ; ig<1 ; ig++ ) {

        gasP_t *p1 = &d->gasP[ig];
        // for testing, just do one octant
        /*if (    p1->r[0]>0. && p1->r[0]<0.01 &&
                p1->r[1]>0. && p1->r[1]<0.01 ) {*/
        
        std::vector< double > p1rhat (3);
        double norm2=0;
        for ( int jj=0 ; jj<3 ; jj++ ) {
            norm2+=square(p1->r[jj]);
        }

        double norm = sqrt(norm2);
    
        for ( int jj=0 ; jj<3 ; jj++ ) {
            p1rhat[jj] = p1->r[jj]/norm;
        }
    
        for ( int jg = 0 ; jg<d->ng ; jg++ ) {
            if ( ig!=jg ) {
                gasP_t *p2 = &d->gasP[jg];
            
                double rad2=0.;
                for ( int jj=0 ; jj<3 ; jj++ ) {
                    rad2+=square(p2->r[jj]);
                }
                //double h2 = square(p2->h);
                //std::cout << rad2 << "," << norm2 << std::endl;
                if ( rad2<norm2 && rad2>1.e-8) {
                    double dotprod=0.;
                    for ( int jj=0 ; jj<3 ; jj++ ) {
                        dotprod+=p2->r[jj]*p1->r[jj];
                    }
            
                    if ( dotprod>0. ) {
                    
                    /*if (    p2->r[0]>0. && p2->r[0]<0.01 &&
                            p2->r[1]>0. && p2->r[1]<0.01 ) {*/
                        double ldist = sqrt(
                            square(p2->r[1]*p1rhat[2] - p1rhat[1]*p2->r[2]) +
                            square(p2->r[2]*p1rhat[0] - p1rhat[2]*p2->r[0]) +
                            square(p2->r[0]*p1rhat[1] - p1rhat[0]*p2->r[1])
                                            );
                        //double dist = fabs(p2->r[0]*p1rhat[1] - p1rhat[0]*p2->r[1]);
                        //std::cout << dist << "-" << p2->h << std::endl;
                        if ( ldist<=p2->h ) {
                            double contrib = kernel->flat_w(ldist,p2->h)*p2->m;
                            dep[ig]+=contrib;
//                             myfile << jg << " " << dep[ig] << " " << contrib << " " << ldist  << " " << p2->h << " " << rad2 << " " << p2->r[0] << " " << p2->r[1] << " " << p2->r[2] << " " << atan2(p2->r[1],p2->r[0]) << std::endl;
                            //dep[ig]+=1;
                        }
                        //dep[ig]+=1;
                    }
                }
            }
        }
        //}
//         myfile.close();
    }
    return dep;
}

std::vector< double > depthrays::SPH_depth2(gizData_t *d) {
    std::vector< double > dep ( d->ng,0.);
//   OLD VERSION: 
//     Bin particles by angle, plus a ``central bin'' for particles whose
//     smoothing lengths are smaller than the distance between angle bins
//     and who need to be checked by *everyone*
//     theta bins are 1-n, ``central'' bin is 0
//   NEW VERSION:
//     Bin particles by angle, and give a ``width'' for particles
//     so we know how many adjacent bins to check

    std::clock_t start;
    double duration;

    start = std::clock();

// TODO - bin by radius too?
    int n_theta_bins = 1024;
    double bin_angular_width = 2.*M_PI/n_theta_bins;
    double bin_angular_width2 = square(bin_angular_width);
    //std::vector< std::vector <int> > theta_bins(n_theta_bins+1, std::vector< int >(d->ng,0) );
    //int nbin[n_theta_bins+1];
    std::vector< std::vector <int> > theta_bins(n_theta_bins, std::vector< int >(d->ng,0) );
    int nbin[n_theta_bins];
    std::ofstream myfile;
    myfile.open("bintest.dat");
    
    std::vector< int > pbins(d->ng,0); // what bin this particle is in
    //std::vector< int > pbin_w(d->ng,0); // how many bins sideways this particle crosses into
    
    //for ( int ii=0 ; ii<n_theta_bins+1 ; ii++ ) {
    for ( int ii=0 ; ii<n_theta_bins ; ii++ ) {
        nbin[ii] = 0;
    }

    for ( int ig=0 ; ig<d->ng ; ig++ ) {
        gasP_t *p1 = &d->gasP[ig];

        double theta = atan2(p1->r[1],p1->r[0]);
        double rad2 = 0;
        for ( int jj=0 ; jj<2 ; jj++ ) {
            rad2+=square(p1->r[jj]);
        }

        double physical_bin_width2 = rad2*bin_angular_width2;
        double h2 = square(p1->h);

        int ibin, bin_w;
        ibin = ((theta/2./M_PI)+.5)*n_theta_bins;
        // deal with rounding errors
        if ( ibin<0 ) ibin=0;
        if ( ibin>=n_theta_bins ) ibin=n_theta_bins;
        
        if ( rad2<=h2 ) {
            bin_w = n_theta_bins/2; // central particles - look at all bins
        } else {
            if ( physical_bin_width2>h2 ) {
                bin_w = 1; // small particles, only look at adjacent bins
            } else { //
                bin_w = ceil((sqrt(h2/physical_bin_width2)));
                //std::cout << "T" << bin_w << " " << h2 << " " << physical_bin_width2 << " " << ceil(2.*sqrt(h2/physical_bin_width2)) << " " << (2.*sqrt(h2/physical_bin_width2))  << std::endl;
                if ( bin_w<1 ) bin_w = 1;
                if ( bin_w>n_theta_bins/2 ) bin_w = n_theta_bins/2;
            }
        }
        
        /*if ( ig==10668 || ig==0 ) {
            std::cout << "X" << ig << " " << bin_w << " " << p1->h << " " << ibin << " " << theta << " " << sqrt(rad2) << " " << bin_angular_width << std::endl;
        }*/
        
        //bin_w = n_theta_bins/2; // for testing, everyone checks everything
        
        /*if ( physical_bin_width2<=h2 ) {
            ibin = 0;
        } else {
            ibin = ((theta/2./M_PI)+.5)*n_theta_bins+1;

            // deal with rounding errors
            if ( ibin<1 ) ibin=1;
            if ( ibin>n_theta_bins ) ibin=n_theta_bins;
        }*/
        myfile << p1->r[0] << " " << p1->r[1] << " " << p1->r[2] << " " << theta << " " << ibin << " " << bin_w << std::endl;
        
        //std::cout << ibin << " " << nbin[ibin] << std::endl;
        
        pbins[ig] = ibin;
        if ( bin_w==n_theta_bins/2 ) {
            for ( int jj=0 ; jj<n_theta_bins ; jj++ ) {
                theta_bins[jj][nbin[jj]]=ig;
                nbin[jj]++;
            }
        } else {
            for ( int jj=-bin_w ; jj<=bin_w ; jj++ ) {
                int jbin = (jj+ibin+n_theta_bins)%n_theta_bins; // positive mod
                theta_bins[jbin][nbin[jbin]]=ig;
                nbin[jbin]++;
            }
        }
    }
    myfile.close();
    int nbinmax = 0;
    //for ( int ii=0 ; ii<n_theta_bins+1 ; ii++ ) {
    for ( int ii=0 ; ii<n_theta_bins ; ii++ ) {
        //std::cout << nbin[ii] << std::endl;
        nbinmax = std::max(nbin[ii],nbinmax);
    }
    std::cout << "Efficiency:" << ((double)d->ng)/(nbinmax*n_theta_bins) << std::endl;
    std::cout << "Speedup:" << ((double)d->ng)/(nbinmax) << std::endl;
    //exit(0);

    /*for ( int ig = 0 ; ig<d->ng ; ig++ ) {
        myfile << pbins[ig] << std::endl;
    }
    myfile.close();*/
    
    // "scatter" form, with bins

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"decomposition: "<< duration <<'\n';

    start = omp_get_wtime();
//#pragma omp parallel for shared(dep,d,pbins,nbin,theta_bins) default(none) schedule(dynamic,100)
    for ( int ig = 0 ; ig<d->ng ; ig++ ) {
    //for ( int ig = 0 ; ig<1000 ; ig++ ) {
    //for ( int ig = 0 ; ig<1 ; ig++ ) {
        //if ( ig%1000==0 ) std::cout << ig << std::endl;
        
        //std::ofstream myfile;
        //myfile.open("contrib2.dat");


        gasP_t *p1 = &d->gasP[ig];
        
        std::vector< double > p1rhat (3);
        double norm2=0;
        for ( int jj=0 ; jj<3 ; jj++ ) {
            norm2+=square(p1->r[jj]);
        }

        double norm = sqrt(norm2);
    
        for ( int jj=0 ; jj<3 ; jj++ ) {
            p1rhat[jj] = p1->r[jj]/norm;
        }
        
        int ibin = pbins[ig];
//         std::cout << "BIN:" << ibin << std::endl;
        
        // now only need to look through one bin
        // because the distribution is done above
        
        // create list of bins to look through
        /*std::vector< int > binlist(1,0);
        binlist[0] = 0; // always do inner stuff
        if ( ibin==0 ) {
            // just do everything for now
            for ( int ii=1 ; ii<=n_theta_bins ; ii++ ) {
                binlist.push_back(ii);
            }
        } else {
            // only do this bin and adjacent bins
            binlist.push_back(ibin-1);
            if ( binlist[1]==0 ) binlist[1] = n_theta_bins;
            binlist.push_back(ibin);
            binlist.push_back(ibin+1);
            if ( binlist[3]==n_theta_bins+1 ) binlist[3] = 1;
        }*/
        /*std::cout << ibin << " " << binlist.size() << std:: endl;
        for ( int ii=0 ; ii<binlist.size() ; ii++ ) {
            std::cout << binlist[ii] << " " ;
        }
        std::cout << std::endl;
        exit(0);*/
        
        //int cc = 0;
        /*for ( int ii = 0 ; ii<binlist.size() ; ii++ ) {
            int jbin = binlist[ii];*/
        
            //for ( int jg = 0 ; jg<d->ng ; jg++ ) {
            int jbin = ibin;
            for ( int jdex = 0 ; jdex<nbin[jbin] ; jdex++ ) {
                int jg = theta_bins[jbin][jdex];
                if ( ig!=jg ) {
                    //cc++;
                    gasP_t *p2 = &d->gasP[jg];
            
                    double rad2=0.;
                    for ( int jj=0 ; jj<3 ; jj++ ) {
                        rad2+=square(p2->r[jj]);
                    }
                    if ( rad2<norm2 && rad2>1.e-8) {
                        double dotprod=0.;
                        for ( int jj=0 ; jj<3 ; jj++ ) {
                            dotprod+=p2->r[jj]*p1->r[jj];
                        }
            
                        if ( dotprod>0. ) {
                    
                            double ldist2 = (
                                square(p2->r[1]*p1rhat[2] - p1rhat[1]*p2->r[2]) +
                                square(p2->r[2]*p1rhat[0] - p1rhat[2]*p2->r[0]) +
                                square(p2->r[0]*p1rhat[1] - p1rhat[0]*p2->r[1])
                                                );
                            double ldist = sqrt(ldist2);
                            //double h2 = square(p2->h);
                            //if ( ldist2<=h2 ) {
                            if ( ldist<=p2->h ) {
                                //dep[ig]+=kernel->flat_w2_q(ldist2,h2)*p2->m;
                                double contrib = kernel->flat_w(ldist,p2->h)*p2->m;
                                dep[ig]+=contrib;
                                // myfile << jg << " " << dep[ig] << " " << contrib << " " << ldist  << " " << p2->h << " " << rad2 << " " << p2->r[0] << " " << p2->r[1] << " " << p2->r[2] << " " << atan2(p2->r[1],p2->r[0]) << std::endl;
                            }
                        }
                    }
                }
            //}
        }
        //std::cout << ig << " " << cc << std::endl;
        //if ( ig==50 ) exit(0);
        // myfile.close();
    }
    duration = omp_get_wtime()-start;
    std::cout<<"main loop: "<< duration <<'\n';
    return dep;
}