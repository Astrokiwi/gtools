fof:
	f2py --opt=-O2 --f90flags="-O2 -Wall -mcmodel=medium" -c fof.f90 -m fof


# FORT = h5fc
# 
# HDF5INCL =  -I/usr/local/hdf5/include
# # -DH5_USE_16_API
# HDF5LIB  =  -L/usr/local/hdf5/lib -lhdf5 -lz
# 
# LIBS   = -lm $(HDF5LIB)
# INCL = $(HDF5INCL)
# 
# OPTIONS =  -O3 -ldl -Wall -mcmodel=medium
# 
# FFLAGS = $(LIBS) $(INCL) $(OPTIONS) 
# 
# .SUFFIXES : .o .f90 .F
# 
# 
# HDF5TEST_OBJ=gizmo_hdf5_input.o gizmo_hdf5_input_test.o
# HDF5TEST_EXEC=gizmo_hdf5_input
# 
# input: ${HDF5TEST_OBJ}
# 	${FORT} ${HDF5TEST_OBJ} ${FFLAGS} -o ${HDF5TEST_EXEC}
# 	#mv ${HDF5TEST_EXEC} ..
# 
# 
# 
# .f90.o :
# 	${FORT} ${FFLAGS} -c $<
# 
# .F.o :
# 	${FORT} ${FFLAGS} -c $<
# 
# 
# clean :
# 	rm *.o *.mod
