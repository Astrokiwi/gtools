#include <vector>
#include "hdf5.h"
#include <stdio.h>
#include <math.h>


int main() {
    int flashFile=0;
    flashFile = H5Fopen("/export/1/djw/gizmo_public/disc_nocool_midres_out/snapshot_096.hdf5",H5F_ACC_RDONLY,H5P_DEFAULT);
    printf("%d\n",flashFile);



	hid_t dataset,dataspace,status,memspace;
	int rank;
	hsize_t dims_out,max_dims;
	
	int data_size;

	dataset = H5Dopen(flashFile, "/PartType0/ParticleChildIDsNumber",H5P_DEFAULT);
	dataspace = H5Dget_space(dataset);
	rank = H5Sget_simple_extent_ndims(dataspace);
	printf("%d\n",rank);
	status = H5Sget_simple_extent_dims(dataspace,&dims_out,&max_dims);
    printf("%d %d\n",dims_out,max_dims);
	
	H5Sclose(dataspace);
	H5Dclose(dataset);

    hsize_t ng = dims_out;
    //std::vector<double> rho(ng);
    double *rho = new double[ng];

	dataset = H5Dopen(flashFile, "/PartType0/Density",H5P_DEFAULT);
	dataspace = H5Dget_space(dataset);
	
    rank = 1;
	memspace = H5Screate_simple(rank,&ng,NULL);
	status = H5Dread(dataset,H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,&rho[0]);
	
	H5Sclose(memspace);
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	for ( int ii=0 ; ii<ng ; ii+=100 ) {
	    printf("%g ",log10(rho[ii]));
	}
	printf("\n");

    return 1;
}
