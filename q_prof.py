print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys

from sys import path
path.append("src/")

import gizmo_tools

def bin_func(val_p,slice,bin_indices,nbins,func):
    bin_func_val = np.zeros((nbins))
    for ibin in range(nbins):
        if np.sum(bin_indices==ibin)>0:
            bin_func_val[ibin] = func((val_p[slice])[bin_indices==ibin])
        else:
            bin_func_val[ibin] = 0.
    return bin_func_val

print("Running")


run_id = sys.argv[1]
output_dir = sys.argv[2]
snap_str = sys.argv[3]

Nbins = 50

gizmoDir = gizmo_tools.getGizmoDir()
fullDir = gizmoDir+"/"+run_id+"/"+output_dir


f = h5py.File(fullDir+"/snapshot_"+snap_str+".hdf5","r")

xyz = np.array(f["/PartType0/Coordinates"])


u_p = np.array(f["/PartType0/InternalEnergy"])
rho_p = np.array(f["/PartType0/Density"])
vel_p = np.array(f["/PartType0/Velocities"])
m_p = np.array(f["/PartType0/Masses"])



xyz*=1.e3 # from kpc to pc
rho_p*=6.77e-22 # to g/cm**3
u_p*=1.e10 # to erg/g
vel_p*=1.e5 # to cm/s
m_p*=1.989e+43 # 10^10 solar masses to g

G = 6.67259e-8 # in cgs
cs_p = np.sqrt(10.*u_p/9.) # cm/s, assuming gamma=5/3

rad_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2)
rad2d_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
z_p = np.abs(xyz[:,2])

vrad_p = (xyz[:,0]*vel_p[:,0]+xyz[:,1]*vel_p[:,1]+xyz[:,2]*vel_p[:,2])/rad_p
vcirc_p = (xyz[:,1]*vel_p[:,0]-xyz[:,0]*vel_p[:,1])/rad2d_p

omega_p = vcirc_p/(rad2d_p*3.085678e18)

N_part = rad_p.size


print(np.min(u_p))
print(np.min(cs_p))

r_cut = 50. # in pc
vrad_cut = 50. # in km/s
low_u = 4.e9 # about 30 K

disc_slice = (vrad_p<vrad_cut) & (rad_p<r_cut)
ndisc = np.sum(disc_slice)

rad2d_p = rad2d_p[disc_slice]

nbins = 100
rbins = np.sort(rad2d_p)[::ndisc//(nbins+1)]

bin_indices = np.digitize(rad2d_p,rbins)-1

cs_bin = bin_func(cs_p,disc_slice,bin_indices,nbins,np.mean)
omega_bin = bin_func(omega_p,disc_slice,bin_indices,nbins,np.mean)
vcirc_bin = bin_func(vcirc_p,disc_slice,bin_indices,nbins,np.mean)
rad2d_bin = bin_func(rad2d_p,[True]*ndisc,bin_indices,nbins,np.mean)
mass_bin = bin_func(m_p,disc_slice,bin_indices,nbins,np.sum)
area_bin = np.pi*(rbins[1:nbins+1]**2-rbins[:nbins]**2)*3.085678e18**2
surface_density_bin = mass_bin/area_bin

omega2_bin = omega_bin**2
omega2diff_bin = np.gradient(omega2_bin)/np.gradient(rad2d_bin)

kappa_bin = np.sqrt(rad2d_bin*omega2diff_bin+4.*omega2_bin)


Q_bin = cs_bin*kappa_bin/(np.pi*G*surface_density_bin)
Q_lowtemp = np.sqrt(10./9.*low_u)*kappa_bin/(np.pi*G*surface_density_bin)

output_data = [rad2d_bin,Q_bin,Q_lowtemp,cs_bin,kappa_bin,surface_density_bin]
output_data = np.array(output_data).T

np.savetxt("data/toomreqprof"+run_id+output_dir+snap_str+".dat",output_data)






