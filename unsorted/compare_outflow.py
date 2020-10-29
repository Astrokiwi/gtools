import sys
import gizmo_tools
import numpy as np
import h5py
from scipy import interpolate
from scipy import signal
from multiprocessing import Pool
import itertools as it
import sys

G_kms_pc_msun = 0.0043022682

def v_esc(r,m_bh,m_hern,a_hern):
    return np.sqrt(2*G_kms_pc_msun*(m_bh/r + m_hern/(a_hern+r)))


def find_outflow_p(run_id,output_dir,snap_str):
    gizmoDir = gizmo_tools.getGizmoDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    f = h5py.File(fullDir+"/snapshot_"+snap_str+".hdf5","r")

    xyz_p = np.array(f["/PartType0/Coordinates"])
    vel_p = np.array(f["/PartType0/Velocities"])
    id_p = np.array(f["/PartType0/ParticleIDs"])

    xyz_p*=1.e3 # to pc
    rad3d_p = np.sqrt(np.sum(xyz_p**2,axis=1))

    vmag = np.sqrt(np.sum(vel_p**2,axis=1))
    
    outflowing = (v_esc(rad3d_p,1.e6,1.e9,250.)<vmag)
    
    return id_p[outflowing],xyz_p[outflowing,:],rad3d_p[outflowing],vmag[outflowing]


if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    snap0_str = sys.argv[3]
    snap1_str = sys.argv[4]
    
    outflow0_ids,outflow0_locs,rad0,vmag0 = find_outflow_p(run_id,output_dir,snap0_str)
    outflow1_ids,outflow1_locs,rad1,vmag1 = find_outflow_p(run_id,output_dir,snap1_str)
    
    np.savetxt("data/outflow0_locs.dat",np.vstack([outflow0_locs.T,rad0,vmag0]).T)
    np.savetxt("data/outflow1_locs.dat",np.vstack([outflow1_locs.T,rad1,vmag1]).T)
    
    outflow_lost_slice = ~np.isin(outflow0_ids,outflow1_ids)
    outflow_gained_slice = ~np.isin(outflow1_ids,outflow0_ids)
    
    print(np.sum(outflow_lost_slice),np.sum(outflow_gained_slice))
    
    np.savetxt("data/outflow_lost.dat",np.vstack([outflow0_locs[outflow_lost_slice].T,rad0[outflow_lost_slice],vmag0[outflow_lost_slice]]).T)
    np.savetxt("data/outflow_gained.dat",np.vstack([outflow1_locs[outflow_gained_slice].T,rad1[outflow_gained_slice],vmag1[outflow_gained_slice]]).T)
    

