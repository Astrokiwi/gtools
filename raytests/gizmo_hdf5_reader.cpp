// currently only reads gas

#include "hdf5.h"
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>

#include "prototypes.h"


gizmo_hdf5_reader::gizmo_hdf5_reader(char const *inName) {
	gizFile = H5Fopen(inName,H5F_ACC_RDONLY,H5P_DEFAULT);
	
	if( (int) gizFile < 0)
	{
		printf("ERROR: FILE NOT FOUND \n");
		exit(0);
	}
	
	readDims();

	if( (int) gizFile < 0)
	{
		printf("ERROR: FILE HAS INCORRECT DIMENSION \n");
		exit(0);
	}
}

void gizmo_hdf5_reader::readGas() {
    //std::clock_t start;
//     double duration;

    //start = std::clock();
    if ( false )
    {
        double* data;
    
        data = readDoubleArray("/PartType0/Density");
    
        for ( int ig=0 ; ig<myData.ng ; ig++ ) {
            myData.gasP[ig].rho = data[ig];
        }
    
    
        data = readDoubleArray("/PartType0/InternalEnergy");
    
        for ( int ig=0 ; ig<myData.ng ; ig++ ) {
            myData.gasP[ig].uint = data[ig];
        }

        data = readDoubleArray("/PartType0/Masses");
    
        for ( int ig=0 ; ig<myData.ng ; ig++ ) {
            myData.gasP[ig].m = data[ig];
        }

    
        data = readDoubleArray("/PartType0/SmoothingLength");
    
        for ( int ig=0 ; ig<myData.ng ; ig++ ) {
            myData.gasP[ig].h = data[ig];
        }

        data = readDoubleArray("/PartType0/Velocities");
    
        for ( int ig=0 ; ig<myData.ng ; ig++ ) {
            for ( int jj=0 ; jj<3 ; jj++ ) {
                myData.gasP[ig].v[jj] = data[ig*3+jj];
            }
        }

        data = readDoubleArray("/PartType0/Coordinates");
    
        for ( int ig=0 ; ig<myData.ng ; ig++ ) {
            for ( int jj=0 ; jj<3 ; jj++ ) {
                myData.gasP[ig].r[jj] = data[ig*3+jj];
            }
        }
        
        delete data;
    
    }
    
    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //std::cout<<"one array  : "<< duration <<'\n';
    
    
    //start = std::clock();
    {
        double *data[6];
    
        data[0] = readDoubleArray("/PartType0/Density");
        data[1] = readDoubleArray("/PartType0/InternalEnergy");
        data[2] = readDoubleArray("/PartType0/Masses");
        data[3] = readDoubleArray("/PartType0/SmoothingLength");
        data[4] = readDoubleArray("/PartType0/Velocities");
        data[5] = readDoubleArray("/PartType0/Coordinates");
    
        for ( int ig=0 ; ig<myData.ng ; ig++ ) {
            myData.gasP[ig].rho = data[0][ig];
            myData.gasP[ig].uint = data[1][ig];
            myData.gasP[ig].m = data[2][ig];
            myData.gasP[ig].h = data[3][ig];
            for ( int jj=0 ; jj<3 ; jj++ ) {
                myData.gasP[ig].v[jj] = data[4][ig*3+jj];
                myData.gasP[ig].r[jj] = data[5][ig*3+jj];
            }
        }
        
        for ( int jj=0 ; jj<6 ; jj++ ) {
            //std::cout << sizeof(data[jj]) << " " <<sizeof(*(data[jj])) << std::endl;
            delete data[jj];
        }
    }
    
    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //std::cout<<"many arrays: "<< duration <<'\n';
}

void gizmo_hdf5_reader::readDims() {
    hid_t dataset,dataspace;//,status;
	//int rank;
	hsize_t dims_out,max_dims;

	dataset = H5Dopen(gizFile, "/PartType0/ParticleChildIDsNumber",H5P_DEFAULT);
	dataspace = H5Dget_space(dataset);
	H5Sget_simple_extent_ndims(dataspace);
	//printf("%d\n",rank);
	H5Sget_simple_extent_dims(dataspace,&dims_out,&max_dims);
    //printf("%d %d\n",dims_out,max_dims);
	
	myData.ng = dims_out;
	
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	myData.gasP.resize(myData.ng);
}

// double* gizmo_hdf5_reader::readDoubleArray(std::string *data_name) {
//     return this->readDoubleArray(data_name->c_str());
// }
 
double* gizmo_hdf5_reader::readDoubleArray(char const *data_name) {
//double* gizmo_hdf5_reader::readDoubleArray(std::string *data_name_in) {
    // copied from flash_reader
	hid_t dataset,dataspace,memspace;//,status;
	//char* data_name = data_name_in->c_str();
	int rank;
	//hsize_t *dims_out,*max_dims;

	double* data;
	
	int data_size;

	if( (int) gizFile < 0)
	{
		printf("CAN NOT READ %s, FILE NOT OPEN\n",data_name);
		exit(0);
	}

	dataset = H5Dopen(gizFile, data_name,H5P_DEFAULT); // Works on others
	dataspace = H5Dget_space(dataset);
	rank = H5Sget_simple_extent_ndims(dataspace);
	//dims_out = new hsize_t[rank];
	//max_dims = new hsize_t[rank];
	
	hsize_t dims_out[rank], max_dims[rank];
	
	H5Sget_simple_extent_dims(dataspace,dims_out,max_dims);
	
	data_size = 1;
	for ( int i=0; i<rank ; i++ )
	{
		data_size*=(int)dims_out[i];
	}
	
//	data = (double *) (malloc(data_size*sizeof(double)));
	data = new double[data_size];
	
	memspace = H5Screate_simple(rank,dims_out,NULL);
	H5Dread(dataset,H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,data);
	H5Sclose(memspace);
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	//delete max_dims;
	//delete dims_out;
	
	return data;

    /*int rank = 1;
    hid_t dataset,dataspace,status,memspace;
    hsize_t memsize;
    
    double *data = new double[myData.ng];

	dataset = H5Dopen(gizFile, data_name,H5P_DEFAULT);
	dataspace = H5Dget_space(dataset);
	
	memsize = (hsize_t)myData.ng;
	
	memspace = H5Screate_simple(rank,&memsize,NULL);
	status = H5Dread(dataset,H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,&data[0]);
	
	H5Sclose(memspace);
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	return data;*/
}

/*double* gizmo_hdf5_reader::readDoubleArray3D(char *data_name) {
    int rank;
    hid_t dataset,dataspace,status,memspace;
    hsize_t memsize;
    
    double *data = new double[myData.ng,3];

	dataset = H5Dopen(gizFile, data_name,H5P_DEFAULT);
	dataspace = H5Dget_space(dataset);
	rank = H5Sget_simple_extent_ndims(dataspace);

	memsize = (hsize_t)(myData.ng);
	
	memspace = H5Screate_simple(rank,&memsize,NULL);
	status = H5Dread(dataset,H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,&data[0]);
	
	H5Sclose(memspace);
	H5Sclose(dataspace);
	H5Dclose(dataset);
	
	return data;
}*/


gizData_t* gizmo_hdf5_reader::getData() {
    return &myData;
}
