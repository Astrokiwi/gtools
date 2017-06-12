print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys

from sys import path
path.append("src/")
import tab_interp

print("Loading table (short form)")
chTab = tab_interp.CoolHeatTab(("coolheat_tab_marta/shrunk_table_labels_080517.dat"),("coolheat_tab_marta/shrunk_table_080517.dat"))
interpTabVec = np.vectorize(chTab.interpTab)

print("Running")


run_id = sys.argv[1]
output_dir = sys.argv[2]
snap_str = sys.argv[3]

Nbins = 50

f = h5py.File("/export/1/djw/gizmos/"+run_id+"/"+output_dir+"/snapshot_"+snap_str+".hdf5","r")

xyz = np.array(f["/PartType0/Coordinates"])
rad_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2)
rad2d_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
z_p = np.abs(xyz[:,2])

N_part = rad_p.size

u_p = np.array(f["/PartType0/InternalEnergy"])
rho_p = np.array(f["/PartType0/Density"])
vel_p = np.array(f["/PartType0/Velocities"])
a_p = np.array(f["/PartType0/Acceleration"])
h_p = np.array(f["/PartType0/SmoothingLength"])
m_p = np.array(f["/PartType0/Masses"])
dt_p = np.array(f["/PartType0/TimeStep"])


rho_p*=6.77e-22 # to g/cm**3
u_p*=1.e10 # to erg/g
rad_p*=1. # already in kpc
vel_p*=1. # in km/s
h_p*=1. # already in kpc
m_p*=1.989e+43 # 10^10 solar masses to g
radaccel_p*=3.0857e21/3.08568e+16**2 # to cm/s/s
a_p*=3.0857e21/3.08568e+16**2 # to cm/s/s
#dt_p*=3.08568e+16 # to seconds

G = 6.67259e-8 # in cgs
cs_p = np.sqrt(u_p)
#a_p*=# 

vrad_p = (xyz[:,0]*vel_p[:,0]+xyz[:,1]*vel_p[:,1]+xyz[:,2]*vel_p[:,2])/rad_p
vcirc_p = (xyz[:,1]*vel_p[:,0]-xyz[:,0]*vel_p[:,1])/rad2d_p

# in km/s/Gyr, roughly
arad_p = (xyz[:,0]*a_p[:,0]+xyz[:,1]*a_p[:,1]+xyz[:,2]*a_p[:,2])/rad_p
radrad_p = (xyz[:,0]*radaccel_p[:,0]+xyz[:,1]*radaccel_p[:,1]+xyz[:,2]*radaccel_p[:,2])/rad_p
acirc_p = (xyz[:,1]*a_p[:,0]-xyz[:,0]*a_p[:,1])/rad2d_p

lambda_cool = 1.e-24 # for testing

molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16

depth_p/=(molecular_mass*proton_mass_cgs) # N in cm**(-2)

mJ_p = np.pi**(2.5)*cs_p**3*G**(-1.5)*rho_p**(-.5)/6.

nH_p = rho_p/(molecular_mass*proton_mass_cgs)

TK_p = gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*u_p
