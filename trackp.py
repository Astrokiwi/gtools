print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
import h5py

import re

print("Running")

run_id = sys.argv[1]
output_dir = sys.argv[2]
i_track = int(sys.argv[3])

if ( len(sys.argv)>4 ):
    snapf = int(sys.argv[4])
else:
    snapf = 0

snapi = 0
fullDir = "/export/1/djw/gizmos/"+run_id+"/"+output_dir


#determine which files to look at

def sort_nicely( l ):
    """ Sort the given list in the way that humans expect.
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )

if ( snapf==0):
    fnames = os.listdir(fullDir)
    sort_nicely(fnames)
    fnames = np.array(fnames)
    #fnames.sort()
    snapshotfilebools = np.array([x.startswith("snapshot") for x in fnames])
    snapshotfiles = fnames[snapshotfilebools]

    ctime = os.path.getmtime(fullDir+"/snapshot_000.hdf5")
    for fname in snapshotfiles[1:]:
        new_snapf = int(fname[9:len(fname)-5])
        new_ctime = os.path.getmtime(fullDir+"/"+fname)
        if ( new_ctime>ctime ) :
            ctime = new_ctime
            snapf = new_snapf

all_outp = []

for snapx in range(snapi,snapf+1):
    infile = fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5"
    f = h5py.File(infile,"r")
    
    header = f["/Header"]
    time = header.attrs.get("Time")
    time*= 0.9778e9 # to yr
    time/=1.e6 # to Myr

    xyz = np.array(f["/PartType0/Coordinates"])
    rad_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2)
    rad2d_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
    z_p = np.abs(xyz[:,2])


    u_p = np.array(f["/PartType0/InternalEnergy"])
    rho_p = np.array(f["/PartType0/Density"])
    vel_p = np.array(f["/PartType0/Velocities"])
    a_p = np.array(f["/PartType0/Acceleration"])
    h_p = np.array(f["/PartType0/SmoothingLength"])
    m_p = np.array(f["/PartType0/Masses"])
    dt_p = np.array(f["/PartType0/TimeStep"])

    if "/PartType0/AGNHeat" in f:
        agn_heat_p = np.array(f["/PartType0/AGNHeat"])
    else:
        agn_heat_p = np.zeros(N_part)

    if "/PartType0/RadiativeAcceleration" in f:
        radaccel_p = np.array(f["/PartType0/RadiativeAcceleration"])
    else:
        radaccel_p = np.zeros((N_part,3))

    if "/PartType0/AGNColDens" in f:
        depth_p = np.array(f["/PartType0/AGNColDens"]) # surface density units
    else:
        depth_p = np.zeros(N_part)

    agn_heat_p*=1e10/3.08568e+16# to erg/s/g

    depth_p*=(1.989e+43/3.086e+21**2) # to g/cm**2

    rho_p*=6.77e-22 # to g/cm**3
    u_p*=1.e10 # to erg/g
    rad_p*=1. # already in kpc
    vel_p*=1. # in km/s
    h_p*=1. # already in kpc
    m_p*=1.989e+43 # 10^10 solar masses to g
    radaccel_p*=3.0857e21/3.08568e+16**2 # to cm/s/s
    a_p*=3.0857e21/3.08568e+16**2 # to cm/s/s

    cs_p = np.sqrt(u_p)

    vrad_p = (xyz[:,0]*vel_p[:,0]+xyz[:,1]*vel_p[:,1]+xyz[:,2]*vel_p[:,2])/rad_p
    vcirc_p = (xyz[:,1]*vel_p[:,0]-xyz[:,0]*vel_p[:,1])/rad2d_p

    # in km/s/Gyr, roughly
    arad_p = (xyz[:,0]*a_p[:,0]+xyz[:,1]*a_p[:,1]+xyz[:,2]*a_p[:,2])/rad_p
    radrad_p = (xyz[:,0]*radaccel_p[:,0]+xyz[:,1]*radaccel_p[:,1]+xyz[:,2]*radaccel_p[:,2])/rad_p
    acirc_p = (xyz[:,1]*a_p[:,0]-xyz[:,0]*a_p[:,1])/rad2d_p

    molecular_mass = 4./(1.+3.*.76)
    proton_mass_cgs = 1.6726e-24
    gamma_minus_one = 5./3.-1.
    boltzmann_cgs = 1.38066e-16

    depth_p/=(molecular_mass*proton_mass_cgs) # N in cm**(-2)

    nH_p = rho_p/(molecular_mass*proton_mass_cgs)
    TK_p = gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*u_p

    id_p = np.array(f["/PartType0/ParticleIDs"])

    tracked_vals = [rad_p,vrad_p,radrad_p,TK_p,agn_heat_p]
    this_p = np.argwhere(id_p==i_track)[0][0]
    
    outp = [time]
    for val_array in tracked_vals:
        outp+=[val_array[this_p]]
    
    all_outp+=[outp]

all_outp = np.array(all_outp)

np.savetxt("data/trackp"+run_id+output_dir+"_"+str(i_track)+".dat",all_outp)





