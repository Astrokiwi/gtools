#include "structs.h"
#include "hdf5.h"
#include <math.h>
#include <list>
#include <unordered_set>
#include <set>
#include <vector>
#include <deque>

class gizmo_hdf5_reader {
    public:
        gizmo_hdf5_reader(char *inName);

        void    readGas();
        
        gizData_t* getData();

    private:
        void    readDims();
		hid_t	gizFile;
		double* readDoubleArray(char *data_name);
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


class rebeamrays {
    public:
        rebeamrays();
        
        std::vector< double > SPH_fluxes_thin(gizData_t *d, std::vector< double >* lums);
        std::vector< double > SPH_fluxes_naive(gizData_t *d, std::vector< double >* lums);
        std::vector< double > SPH_fluxes_grid(gizData_t *d, std::vector< double >* lums);
        std::vector< double > SPH_fluxes_tree(gizData_t *d, std::vector< double >* lums);
        
    private:
        Kernel* kernel;
        
        bool between_particles(gasP_t *p1,gasP_t *p2,gasP_t *p3);
        double intersect_d2_nonorm(gasP_t *p1,gasP_t *p2,gasP_t *p3);
        double get_d12_norm(gasP_t *p1,gasP_t *p2);
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
        
        
        SET_TYPE<int >* beam(double r1[3], double r2[3]);

        SET_TYPE<int >* beam_level(int ilevel, double r1[3], double r2[3]);
        octtreelist* getCellAt(int ilevel,int ix,int iy,int iz);
        PLIST_TYPE<int>* getPList();
        
        void dump();
        void matrix_dump();
        void sort();
        
    private:
        int itreer(double r, int idim, int ilevel);
        double rtreei(int ii, int idim, int ilevel);
        double r_left[3];
        int ixyz[3];
        double size;
        int level;
        octtreelist **nodes;
        PLIST_TYPE<int> *myP;
        int maxlevel;
};

//static inline double square(double x);

static inline double square(double x) {
    return x*x;
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
