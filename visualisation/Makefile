pyread:
	# production
	f2py --opt=-O3 --f90flags="-Wall -mcmodel=large -fopenmp"  --f77flags="-Wall -mcmodel=large -fopenmp" -lgomp -c sph_plotter.f90 -m sph_plotter
	# debug
# 	f2py --opt=-O0 --f90flags="-g -Wall -mcmodel=medium -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow -fopenmp"  --f77flags="-g -Wall -mcmodel=medium -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow -fopenmp" -lgomp -c sph_plotter.f90 -m sph_plotter

test:
    # profiling
# 	gfortran -o ray_test sph_plotter.f90 -Wall -mcmodel=medium -fopenmp -lgomp -pg -g  
    # debuggering
# 	gfortran -O0 -g -o ray_test sph_plotter.f90 -Wall -mcmodel=medium -fopenmp
    # full speed test
	gfortran -O3 -o ray_test sph_plotter.f90 -Wall -mcmodel=medium -fopenmp -lgomp

clean:
	rm *.o *.mod