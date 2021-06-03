# from distutils.core import setup, Extension
import numpy

from numpy.distutils.core import setup, Extension

extension_mod = Extension("gtools._tab_interp",
                        ["gtools/_tab_interp_module.cc", "gtools/dust_temp_interp.cpp"])

extension_f = Extension(name='gtools.sph_plotter',
                 sources=['gtools/sph_plotter.f90'],
                 f2py_options=['--quiet'],
                 extra_f77_compile_args=["-Wall","-mcmodel=medium","-fopenmp","-fallow-invalid-boz","-lgomp"],
                 extra_f90_compile_args=["-Wall","-mcmodel=medium","-fopenmp","-fallow-invalid-boz","-lgomp"]
                )
# 
# 	python -m numpy.f2py --opt=-O3 --f90flags="-Wall -mcmodel=medium -fopenmp -fallow-invalid-boz"  --f77flags="-Wall -mcmodel=medium -fopenmp -fallow-invalid-boz" -lgomp -c sph_plotter.f90 -m sph_plotter
# OPTIONS = -O3 -ldl -Wall -mcmodel=medium


setup(name='gtools',
      version='0.1',
      description='Some plotting tools for gizmo SPH',
#      url='http://github.com/',
      author='David Williamson',
      author_email='d.j.williamson@soton.ac.uk',
      license='GPLv3',
      packages=['gtools'],
      ext_modules=[extension_mod,extension_f],
      include_dirs = [numpy.get_include(),'.','./gtools'])