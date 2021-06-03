all:
	make sph_plotter
	make interplib

sph_plotter:
	# production
	python -m numpy.f2py --opt=-O3 --f90flags="-Wall -mcmodel=medium -fopenmp -fallow-invalid-boz"  --f77flags="-Wall -mcmodel=medium -fopenmp -fallow-invalid-boz" -lgomp -c sph_plotter.f90 -m sph_plotter
	# debug
# 	python -m numpy.f2py --opt=-O0 --f90flags="-g -Wall -mcmodel=medium -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow -fopenmp"  --f77flags="-g -Wall -mcmodel=medium -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow -fopenmp" -lgomp -c sph_plotter.f90 -m sph_plotter



CC       =  gcc
CXX      =  g++

OPTIONS = -O3 -ldl -Wall -mcmodel=medium
#OPTIONS = -O0 -g3 -ldl -Wall -mcmodel=medium

.SUFFIXES : .o .c .cpp .h

.cpp.o :
	${CXX} ${OPTIONS} ${INCL} -c  ${LIBS} $<

.c.o :
	${CC} ${OPTIONS} ${INCL} -c  ${LIBS} $<

interplib:
	swig -python -c++ -o _tab_interp_module.cc tab_interp.i
	#python3 setup.py build_ext --inplace


clean:
	rm *.o *.mod
