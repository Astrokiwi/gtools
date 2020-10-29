print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys

from sys import path
path.append("../src/")

import gizmo_tools

import pynbody
from pynbody.analysis import profile

def bin_func(val_p,slice,bin_indices,nbins,func):
    bin_func_val = np.zeros((nbins))
    for ibin in range(nbins):
        if np.sum(bin_indices==ibin)>0:
            bin_func_val[ibin] = func((val_p[slice])[bin_indices==ibin])
        else:
            bin_func_val[ibin] = 0.
    return bin_func_val

if __name__=='__main__':
    print("Running")


    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    snap_str = sys.argv[3]

    header,snap = gizmo_tools.load_gizmo_nbody("devel_inflow","fullrun","100",load_vals=["temp"])

    p = profile.Profile(snap,min='.0001 pc', max='1 pc')

    omega_bin = p['vtheta']/p['rbins']
    omega2_bin = omega_bin**2
    omega2diff_bin = np.gradient(omega2_bin)/np.gradient(p['rbins'])
    kappa_bin = (p['rbins']*omega2diff_bin+4.*omega2_bin)**(1,2)

#     r_cut = 70. # in pc
#     vrad_cut = 50. # in km/s
#     vrad_cut*=1.e5 # to cm/s
#     low_u = 4.e9 # about 30 K
# 
#     disc_slice = (vrad_p<vrad_cut) & (rad_p<r_cut)
#     ndisc = np.sum(disc_slice)
#     rad2d_p = rad2d_p[disc_slice]

    Q_bin = p['cs']*kappa_bin/(np.pi*pynbody.units.Unit('G')*p['density'])
    print(p['rbins']*omega2diff_bin+4.*omega2_bin)
    print(omega_bin)
    print(omega2_bin)
    print(kappa_bin)
    
    print(Q_bin)

#     output_data = [rad2d_bin,Q_bin,Q_lowtemp,cs_bin,kappa_bin,surface_density_bin]
#     output_data = np.array(output_data).T
# 
#     np.savetxt("data/toomreqprof"+run_id+output_dir+snap_str+".dat",output_data)






