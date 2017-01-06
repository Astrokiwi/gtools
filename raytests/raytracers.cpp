#include "prototypes.h"
#include <vector>
#include <math.h>
#include <ctime>
#include <iostream>

raytracers::raytracers() {
    // default values
    L = 128;
    w = 2.e-2; // in kpc, so 20 pc
    updateSurfaceArea();
    
    kernel = new Kernel();
}
        
std::vector< std::vector<double> > raytracers::SPH_surf(gizData_t *d) {
    std::vector< std::vector<double> > g ( L, std::vector<double> (L,0) );
    
    std::clock_t start;
    double duration;

    start = std::clock();
    
    {
        for ( int ig = 0 ; ig<d->ng ; ig++ ) {
    //    for ( int ig = 500 ; ig<501 ; ig++ ) {
            gasP_t *myP = &d->gasP[ig];
        
            int ih = myP->h/cell_w;
        
            if ( ih<=1 ) {
        
                std::vector<int> ixyz = posToGrid(myP->r);
                if ( inGrid(ixyz) ) {
                    g[ixyz[0]][ixyz[1]]+=myP->m;
                }
            } else {
                std::vector<int> ixyz = posToGrid(myP->r);
                std::vector<int> ixyz_low(2),ixyz_high(2);
            
                for ( int jj=0; jj<2 ; jj++ ) {
                    ixyz_low[jj] = std::max(0,ixyz[jj]-ih);
                    ixyz_high[jj] = std::min(L-1,ixyz[jj]+ih);
                }
            
                double r2;
                for ( int ix=ixyz_low[0] ; ix<=ixyz_high[0] ; ix++ ) {
                    for ( int iy=ixyz_low[1] ; iy<=ixyz_high[1] ; iy++ ) {
                        if ( inGrid(ix,iy) ) {
                            r2=sqrt(pow((ix-ixyz[0])*cell_w,2)+pow((iy-ixyz[1])*cell_w,2));
                            if ( r2<=myP->h ) {
                                g[ix][iy]+=myP->m*kernel->flat_w(r2,myP->h)*cellSurf;
                            }
                        }
                    }
                }
            }
        
        
        }
    }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"sqrt version: "<< duration <<'\n';
    convertSumToColumnDensity(&g);
    
    return g;
}

        
std::vector< std::vector<double> > raytracers::SPH_surf2(gizData_t *d) {
    std::vector< std::vector<double> > g ( L, std::vector<double> (L,0) );
    
    std::clock_t start;
    double duration;

    start = std::clock();
    
    {
        for ( int ig = 0 ; ig<d->ng ; ig++ ) {
    //    for ( int ig = 500 ; ig<501 ; ig++ ) {
            gasP_t *myP = &d->gasP[ig];
        
            int ih = myP->h/cell_w;
        
            if ( ih<=1 ) {
        
                std::vector<int> ixyz = posToGrid(myP->r);
                if ( inGrid(ixyz) ) {
                    g[ixyz[0]][ixyz[1]]+=myP->m;
                }
            } else {
                std::vector<int> ixyz = posToGrid(myP->r);
                std::vector<int> ixyz_low(2),ixyz_high(2);
            
                for ( int jj=0; jj<2 ; jj++ ) {
                    ixyz_low[jj] = std::max(0,ixyz[jj]-ih);
                    ixyz_high[jj] = std::min(L-1,ixyz[jj]+ih);
                }
            
                double r2;
                double h2 = square(myP->h);
                for ( int ix=ixyz_low[0] ; ix<=ixyz_high[0] ; ix++ ) {
                    for ( int iy=ixyz_low[1] ; iy<=ixyz_high[1] ; iy++ ) {
                        if ( inGrid(ix,iy) ) {
                            r2=square((ix-ixyz[0])*cell_w)+square((iy-ixyz[1])*cell_w);
                    
                            if ( r2<=h2 ) {
                                g[ix][iy]+=myP->m*kernel->flat_w2(r2,h2)*cellSurf;
                            }
                        }
                    }
                }
            }
        
        
        }
    }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"square version: "<< duration <<'\n';
    convertSumToColumnDensity(&g);
    
    return g;
}

std::vector< std::vector<double> > raytracers::SPH_surf3(gizData_t *d) {
    std::vector< std::vector<double> > g ( L, std::vector<double> (L,0) );
    
    std::clock_t start;
    double duration;

    start = std::clock();
    
    {
        double cell_w_sq = square(cell_w);
        for ( int ig = 0 ; ig<d->ng ; ig++ ) {
    //    for ( int ig = 500 ; ig<501 ; ig++ ) {
            gasP_t *myP = &d->gasP[ig];
        
            int ih = myP->h/cell_w;
        
            if ( ih<=1 ) {
        
                std::vector<int> ixyz = posToGrid(myP->r);
                if ( inGrid(ixyz) ) {
                    g[ixyz[0]][ixyz[1]]+=myP->m;
                }
            } else {
                //std::vector<int> ixyz = posToGrid(myP->r);
                //std::vector<int> ixyz_low(2),ixyz_high(2);
                int ixyz_low[2],ixyz_high[2];
                int ixyz[2];
            
                for ( int jj=0; jj<2 ; jj++ ) {
                    ixyz[jj] = posToGrid(myP->r[jj]);
                    ixyz_low[jj] = std::max(0,ixyz[jj]-ih);
                    ixyz_high[jj] = std::min(L-1,ixyz[jj]+ih);
                }
            
                double r2;
                double h2 = square(myP->h);
                double invh2 = 1./h2;
                for ( int ix=ixyz_low[0] ; ix<=ixyz_high[0] ; ix++ ) {
                    for ( int iy=ixyz_low[1] ; iy<=ixyz_high[1] ; iy++ ) {
                        if ( inGrid(ix,iy) ) {
                            r2=(square((ix-ixyz[0]))+square((iy-ixyz[1])))*cell_w_sq;
                            
                            if ( r2<=h2 ) {
                                g[ix][iy]+=myP->m*kernel->flat_w2_q(r2,invh2);//*cellSurf;
                            }
                        }
                    }
                }
            }
        
        
        }
    }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"squared circle version: "<< duration <<'\n';
    //convertSumToColumnDensity(&g);
    
    return g;
}


std::vector< std::vector<double> > raytracers::naive_surf(gizData_t *d) {
    std::vector< std::vector<double> > g ( L, std::vector<double> (L,0) );
    
    for ( int ig = 0 ; ig<d->ng ; ig++ ) {
//    for ( int ig = 500 ; ig<501 ; ig++ ) {
        gasP_t *myP = &d->gasP[ig];
        
        std::vector<int> ixyz = posToGrid(myP->r);
        if ( inGrid(ixyz) ) {
            g[ixyz[0]][ixyz[1]]+=myP->m;
        }
    }
    
    convertSumToColumnDensity(&g);
    
    
    
    return g;
}

void raytracers::setL(int L) {
    this->L = L;
    updateSurfaceArea();
}

void raytracers::setw(double w) {
    this->w = w;
    updateSurfaceArea();
}

void raytracers::updateSurfaceArea() {
    cell_w = w/L;
    cellSurf = pow(cell_w,2);
}

void raytracers::convertSumToColumnDensity(std::vector< std::vector<double> > *g) {
    for ( int ix=0 ; ix<L ; ix++ ) {
        for ( int iy=0 ; iy<L ; iy++ ) {
            (*g)[ix][iy]/=cellSurf;
        }
    }
}


bool raytracers::inGrid(int ix, int iy) {
    if ( ix<0 || ix>=L ) {
        return false;
    }
    if ( iy<0 || iy>=L ) {
        return false;
    }
    return true;
}

bool raytracers::inGrid(std::vector<int> ixyz) {
    for ( int jj=0 ; jj<3 ; jj++ ) {
        if ( ixyz[jj]<0 || ixyz[jj]>=L ) {
            return false;
        }
    }
    return true;
}

std::vector<int> raytracers::posToGrid(double *r) {
    std::vector<int> ixyz(3);
    for ( int jj=0 ; jj<3 ; jj++ ) {
        ixyz[jj] = (int)((r[jj]/w+.5)*L);
    }
    return ixyz;
}


int raytracers::posToGrid(double r) {
    return (int)((r/w+.5)*L);;
}

