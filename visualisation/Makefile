pyread:
	# production
	f2py --opt=-O2 --f90flags="-Wall -mcmodel=medium"  --f77flags="-Wall -mcmodel=medium" -c sph_plotter.f90 -m sph_plotter
	# debug
	#f2py --opt=-O0 --f90flags="-g -Wall -mcmodel=medium -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow"  --f77flags="-g -Wall -mcmodel=medium -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow" -c sph_plotter.f90 -m sph_plotter

clean:
	rm *.o *.mod