print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys

i_snap = (sys.argv[1])


f = h5py.File("/export/1/djw/gizmo_public/disc_nocool_out/snapshot_000.hdf5","r")

xyz = np.array(f["/PartType0/Coordinates"])
#rad_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2)
rad_p_0 = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
z_p_0 = np.abs(xyz[:,2])

u_p_0 = np.array(f["/PartType0/InternalEnergy"])
rho_p = np.array(f["/PartType0/Density"])

rho_p*=6.77e-22 # to g/cm**3
u_p*=1.e10 # to erg/g
rad_p*=1. # already in kpc

vel_p = np.array(f["/PartType0/Velocities"])
v_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2)

molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16

print("not finished")