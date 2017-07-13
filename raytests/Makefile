# FORT = gfortran
#FORT = h5fc
CXX = g++

HDF5INCL =  -I/usr/local/hdf5/include
# -DH5_USE_16_API
HDF5LIB  =  -L/usr/local/hdf5/lib -lhdf5 -lz

LIBS   = -lm $(HDF5LIB)
INCL = $(HDF5INCL)

OPTIONS =  -std=c++11 -fopenmp -O3 -ldl -Wall -mcmodel=medium
#OPTIONS =  -std=c++11 -fopenmp -O0 -g3 -ldl -Wall -mcmodel=medium
#OPTIONS =  -pg -std=c++11 -fopenmp -ldl -Wall -mcmodel=medium
#OPTIONS =  -g -std=c++11 -fopenmp -O0 -ldl -Wall -mcmodel=medium
#OPTIONS =  -g -std=c++11 -fno-inline -fopenmp -ldl -Wall -mcmodel=medium
#OPTIONS =  -std=c++11 -fopenmp -fprofile-arcs -ftest-coverage -g3 -ldl -Wall -mcmodel=medium
#OPTIONS =  -std=c++11 -fopenmp -fprofile-arcs -ftest-coverage -g -ldl -Wall -mcmodel=medium

# FFLAGS = $(LIBS) $(INCL) $(OPTIONS) 
#CFLAGS = $(LIBS) $(OPTIONS) ${INCL}

.SUFFIXES : .o .cpp .h

.cpp.o :
	${CXX} ${OPTIONS} ${INCL} -c  ${LIBS} $<


REBEAMOBJS=gizmo_hdf5_reader.o mintemptree.o octtreelist.o gridlist.o rebeamrays.o Kernel.o rebeamcomp.o
rebeam: ${REBEAMOBJS}
	${CXX} ${OPTIONS} ${INCL} ${REBEAMOBJS} ${LIBS} -o rebeamcomp


DEPTHOBJS=gizmo_hdf5_reader.o depthrays.o Kernel.o depthraycomp.o
depthcomp: ${DEPTHOBJS}
	${CXX} ${OPTIONS} ${INCL} ${DEPTHOBJS} ${LIBS} -o depthraycomp


COMPOBJS=gizmo_hdf5_reader.o raytracers.o Kernel.o raycomp.o
raycomp: ${COMPOBJS}
	${CXX} ${OPTIONS} ${INCL} ${COMPOBJS} ${LIBS} -o raycomp

#TESTOBJS=hdf5_read_test.o
TESTOBJS=gizmo_hdf5_reader.o classy_test.o
test: ${TESTOBJS}
	${CXX} ${OPTIONS} ${INCL} ${TESTOBJS} ${LIBS} -o testhdf5


# %.o: %.mod

# .SUFFIXES : .o .f90 .F
# 
# 
# HDF5TEST_OBJ=gizmo_hdf5_input.o gizmo_hdf5_input_test.o
# HDF5TEST_EXEC=gizmo_hdf5_input
# 
# input: ${HDF5TEST_OBJ}
# 	${FORT} ${FFLAGS} -o ${HDF5TEST_EXEC} ${HDF5TEST_OBJ}
# 	#mv ${HDF5TEST_EXEC} ..
# 
# 
# 
# .f90.o :
# 	${FORT} ${FFLAGS} -c $<
# 
# .F.o :
# 	${FORT} ${FFLAGS} -c $<


clean :
	rm *.o
# *.mod
