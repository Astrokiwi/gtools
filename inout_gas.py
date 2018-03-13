print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
from joblib import Parallel, delayed
import h5py

from sys import path
path.append("../")
import gizmo_tools

def summarise_file(infile,threshold_rad):
    with h5py.File(infile,"r") as f:
        header = f["/Header"]
        time = header.attrs.get("Time")
        time*= 0.9778e9 # to yr
        time/=1.e6 # to Myr

        xyz = np.array(f["/PartType0/Coordinates"]) # kpc

        m_p = np.array(f["/PartType0/Masses"]) # 10^10 msun

    m_p*=1e+10 # 10^10 solar masses to solar masses
    xyz*=1.e3 # to pc
    
    
    rad2 = np.sum(xyz[:,0:2]**2,axis=1)
    threshold_rad2 = threshold_rad**2
    
    particles_inside = (rad2<threshold_rad2)
    
    mtot = np.sum(m_p)
    
    mass_inside = np.sum(m_p[particles_inside])
    mass_outside = mtot-mass_inside
    
    return time,mass_inside,mass_outside
    

if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    threshold_rad = float(sys.argv[3])
    nprocs = int(sys.argv[4])

    snapi = 0
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir)

    gizmoDir = gizmo_tools.getGizmoDir()
    movieDir = gizmo_tools.getMovieDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    
    infiles = [fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5" for snapx in range(snapi,snapf+1)]
    
    #print(summarise_file(infiles[0],threshold_rad))

    data = Parallel(n_jobs=nprocs)(delayed(summarise_file)(infile,threshold_rad) for infile in infiles)
    data = np.array(data)
    
    np.savetxt("data/inout"+run_id+output_dir+".dat",data)

