#include "structs.h"
#include "hdf5.h"
#include <math.h>
#include <list>
#include <unordered_set>
#include <set>
#include <vector>
#include <deque>
#include <iostream>

class gizmo_hdf5_reader {
    public:
        gizmo_hdf5_reader(char const *inName);

        void    readGas();
        
        gizData_t* getData();

    private:
        void    readDims();
		hid_t	gizFile;
		double* readDoubleArray(char const *data_name);
        //double* readDoubleArray(std::string *data_name_in);
		//double* readDoubleArray3D(char *data_name);
        
        gizData_t myData;
};

class Kernel {
    public:
        Kernel();
        
        double w(double x);
        double w(double r, double h);
        double flat_w(double r, double h);
        double flat_w2(double r2, double h2);
        double flat_w2_q(double r2, double h2);
    
        void set_L_flat(int L_in, int zL_in);
    private:
        int L_flat,zL;
        std::vector<double> flattened_table;
        std::vector<double> flattened_table2;
        
        void integrate_flat();
        
};

class raytracers {
    public:
        raytracers();
        
        std::vector< std::vector<double> > SPH_surf(gizData_t *d);
        std::vector< std::vector<double> > SPH_surf2(gizData_t *d);
        std::vector< std::vector<double> > SPH_surf3(gizData_t *d);
        std::vector< std::vector<double> > naive_surf(gizData_t *d);
        
        void setL(int L);
        void setw(double w);

    private:
        std::vector<int> posToGrid(double *r);
        int posToGrid(double r);
        bool inGrid(std::vector<int> ixyz);
        bool inGrid(int ix, int iy);
        void updateSurfaceArea();
        void convertSumToColumnDensity(std::vector< std::vector<double> > *g);

    
        int L;
        double w;
        double cellSurf;
        double cell_w;
        Kernel* kernel;
};



class depthrays {
    public:
        depthrays();
        
        std::vector< double > SPH_depth(gizData_t *d);
        std::vector< double > SPH_depth2(gizData_t *d);

    private:
        Kernel* kernel;
};


class gridList {
    public:
        gridList(int L);
        std::list<int >* listAt(int ix, int iy, int iz);
        
    private:
        std::vector<std::list<int > > *listVector;
        int L;
};


class octtreelist {
    public:
        octtreelist(double r_left[3], double size);
        octtreelist(int level, int ixyz[3]);
        
        void addPCoord(int ip,int ilevel,int ix,int iy,int iz);
        void head_addP(int ip, double r_p[3], double h);
        
        
        //SET_TYPE<int >* beam(double r1[3], double r2[3]);
        //SET_TYPE<int >* beam_level(int ilevel, double r1[3], double r2[3]);
        void beam_level(int ilevel, double r1[3], double r2[3],std::vector<PLIST_TYPE<int>*> *list_of_lists);
        std::vector<PLIST_TYPE<int>*>* beam(double r1[3], double r2[3]);

        
        octtreelist* getCellAt(int ilevel,int ix,int iy,int iz);
        octtreelist* getCellAt(int ilevel,int ixyz[3]);
        PLIST_TYPE<int>* getPList();
        
        void dump();
        void matrix_dump();
        void sort();
        
    private:
        int itreer(double r, int idim, int ilevel);
        double rtreei(int ii, int idim, int ilevel);
        PLIST_TYPE<int> myP;
        double r_left[3];
        int ixyz[3];
        double size;
        int level;
        octtreelist **nodes;
        //PLIST_TYPE<int> *myP;
        int maxlevel;
};

class mintemptree {
    public:
        mintemptree(double r_cent[3], double size);
        void addp(double r_p[3], int ip, double temp);
        
        void dump(int ilevel);
        void dump();


        std::vector<int>* getThreshList(gizData_t* d, int ip, double mtot);
        void propagateThreshList(gizData_t* d, int ip,std::vector<int>* threshList,int ilevel, double mtot);

    private:
        void placeInChildNode(double r_p[3], int ip, double temp);
        
        double r_cent[3];
        double size;
        
        mintemptree **nodes;
        
        double r_p[3];
        int ip;
        double minTemp;
        
};




class rebeamrays {
    public:
        rebeamrays();
        
        std::vector< double > SPH_fluxes_thin(gizData_t *d, std::vector< double >* lums);
        std::vector< double > SPH_fluxes_naive(gizData_t *d, std::vector< double >* lums);
        std::vector< double > SPH_fluxes_grid(gizData_t *d, std::vector< double >* lums);
        std::vector< double > SPH_fluxes_tree(gizData_t *d, std::vector< double >* lums);
        
    private:
        Kernel* kernel;

        double one_sph_tree_ray(octtreelist *tree, gasP_t *p1, gasP_t *p2,bool alreadyDone[],gizData_t *d, int ig, int jg, double d12_norm);

        
        //bool between_particles(gasP_t *p1,gasP_t *p2,gasP_t *p3);
        //double intersect_d2_nonorm(gasP_t *p1,gasP_t *p2,gasP_t *p3);
        //double get_d12_norm(gasP_t *p1,gasP_t *p2);
};


//static inline double square(double x);

static inline double square(double x) {
    return x*x;
    //return pow(x,2);
}


static inline double powerfour(double x) {
    return x*x*x*x;
    //return pow(x,2);
}


static inline int square(int x) {
    return x*x;
    //return pow(x,2);
}

// cell location to integer
static inline int igridr(double x, double gridR, int L) {
    return floor(((x+gridR)/2./gridR)*L);
}

// "left" side of cell
static inline double rgridi(int ix, double gridR, int L) {
    return ix*((gridR*2.)/L)-gridR;
}

static inline double intersect_d2_nonorm(gasP_t *p1,gasP_t *p2,gasP_t *p3) {
    double dr31[3],dr32[3];
    double crossp;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        dr31[ii] = p3->r[ii]-p1->r[ii];
        dr32[ii] = p3->r[ii]-p2->r[ii];
    }
    
    crossp = square( dr31[1]*dr32[2]-dr31[2]*dr32[1]);
    crossp+= square(-dr31[0]*dr32[2]+dr31[2]*dr32[0]);
    crossp+= square( dr31[0]*dr32[1]-dr31[1]*dr32[0]);
    
    return crossp;
}




static inline bool between_particles(gasP_t *p1,gasP_t *p2,gasP_t *p3) {
    //bool between = false;    
    double z1=0.,z2=0.,z3=0.;
    // try to avoid cache misses
    /*double dd[3];
    for ( int jj=0 ; jj<3 ; jj++ ) {
        dd[jj] = p1->r[jj]-p2->r[jj];
    }
    for ( int jj=0 ; jj<3 ; jj++ ) {
        z1+=p1->r[jj]*dd[jj];
    }
    for ( int jj=0 ; jj<3 ; jj++ ) {
        z2+=p2->r[jj]*dd[jj];
    }
    for ( int jj=0 ; jj<3 ; jj++ ) {
        z3+=p3->r[jj]*dd[jj];
    }*/
    for ( int jj=0 ; jj<3 ; jj++ ) {
        double dd = p1->r[jj]-p2->r[jj];
        z1+=p1->r[jj]*dd;
        z2+=p2->r[jj]*dd;
        z3+=p3->r[jj]*dd;
    }
    if ( z1>z2 ) {
        if ( z1>z3 && z3>z2 ) {
            //between = true;
            return true;
        }
    } else if ( z1<z2 ) {
        if ( z1<z3 && z3<z2 ) {
            //between = true;
            return true;
        }
    } else {
        // ?? should not happen
        std::cout << "z1=z2 ???" << std::endl;
        exit(0);
    }
    //if ( (z1<z2 && z2<z3) || (z1>z2 && z2>z3 ) ) return true;
    //return between;
    return false;
}


static inline double get_d12_norm(gasP_t *p1,gasP_t *p2) {
    double norm=0.;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        norm+=square(p1->r[ii]-p2->r[ii]);
    }
    
    return norm;
}







// units of area/mass - dust mass fraction * mean dust cross section/mean dust mass, in internal units, i.e. (kpc**2)/(1e10 Msun)
// Here I'm using a dust grain radius of 10 µm and an internal dust density of 6 g/cm**3
// and a dust mass fraction of 1%
static const double dust_physics = .01 * 156.7;

// fake values - these are tabulated later
static const double planck_mean_absorb = .7;
static const double planck_mean_emit_absorb = square(planck_mean_absorb);




