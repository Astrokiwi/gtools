import numpy as np
import h5py
import gizmo_tools
from joblib import Parallel, delayed

vz_cut = 1. # km/s
radcut = 80. # pc
nprocs = 128

def get_rad0(infile):
    with h5py.File(infile,"r") as f:
        xyz_p = np.array(f["/PartType0/Coordinates"]) # kpc
        id_p = np.array(f["/PartType0/ParticleIDs"]).astype(int)
    id_p-=1
    rad_p = np.sqrt(xyz_p[:,0]**2+xyz_p[:,1]**2+xyz_p[:,2]**2)*1.e3
    rad_out = np.zeros_like(rad_p)
    rad_out[id_p] = rad_p
    return rad_out

def get_rad_z(infile,rad0):
    with h5py.File(infile,"r") as f:
        xyz_p = np.array(f["/PartType0/Coordinates"]) # kpc
        vel_p = np.array(f["/PartType0/Velocities"]) # in km/s
        id_p = np.array(f["/PartType0/ParticleIDs"]).astype(int)
        header = f["/Header"]
        time = header.attrs.get("Time")
        time*= 0.9778e9 # to yr
        time/=1.e6 # to Myr
    id_p-=1
    rad_p = np.sqrt(xyz_p[:,0]**2+xyz_p[:,1]**2+xyz_p[:,2]**2)*1.e3
    z_p = np.abs(xyz_p[:,2])*1.e3
    vz_p = np.abs(vel_p[:,2])
    rad0_p = rad0[id_p]
    return time,rad0_p,rad_p,z_p,vz_p


def get_mean_drad_lowvz(infile,rad0):
    time,rad0_p,rad_p,z_p,vz_p = get_rad_z(infile,rad0)
    print(time)
    lowvz_selection = (vz_p<vz_cut)
    notsf_selection = (rad_p<radcut)
    combined_selection = lowvz_selection & notsf_selection
    drad_p = rad_p[combined_selection]-rad0_p[combined_selection]
    return time,np.mean(drad_p)

if __name__ == '__main__':
    run_id = "2014"
    
    for output_dir in ["q2edd05redo","q2edd10_aniso1redo","q2edd10_aniso3redo","q2edd10redo","q2edd20redo","q2redo"]:
#     for output_dir in ["q2edd20redo"]:
    #     output_dir = "q2redo"
    #     output_dir = "q2redo"

        snapi = 0
        snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir)

        gizmoDir = gizmo_tools.getGizmoDir()
        movieDir = gizmo_tools.getMovieDir()
        fullDir = gizmoDir+"/"+run_id+"/"+output_dir

        infiles = [fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5" for snapx in range(snapi,snapf+1)]

        rad0 = get_rad0(infiles[0])


        drad_lowz = Parallel(n_jobs=nprocs)(delayed(get_mean_drad_lowvz)(infile,rad0) for infile in infiles)
    
    #     out_data = [get_mean_drad_lowvz(infile,rad0) for infile in infiles]
        drad_lowz = np.array(drad_lowz)
        print(drad_lowz.shape)
        v = np.gradient(drad_lowz[:,1])/np.gradient(drad_lowz[:,0])
        a = np.gradient(v)/np.gradient(drad_lowz[:,0])
        out_data = np.vstack([drad_lowz.T,v,a]).T
        
        np.savetxt("data/dvz_selected{}{}.dat".format(run_id,output_dir),out_data)
    
    
    
    #     rad0_p,rad_p,z_p,vz_p = get_rad_z(infiles[50],rad0)
    #     
    #     drad_p = rad_p-rad0_p
    #     out_data = np.array([drad_p,z_p,vz_p]).T
    #     
    #     np.savetxt("data/drad_z.dat",out_data)