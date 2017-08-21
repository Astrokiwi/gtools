from distutils.core import setup, Extension

extension_mod = Extension("_tab_interp", ["_tab_interp_module.cc", "dust_temp_interp.cpp"])

setup(name = "tab_interp", ext_modules=[extension_mod])
