import numpy as np
import gizmo_tools
from multiprocessing import Pool
import itertools as it
import h5py

cutrad = 70.e-3
#cutrad2 = cutrad**2
cutzed = 0.01e-3

def load_calc_summary(run_id,output_dir,snap_id):
    gizmoDir = gizmo_tools.getGizmoDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    snap_str = "%03d" % snap_id
    with h5py.File(fullDir+"/snapshot_"+snap_str+".hdf5","r") as f:
        header = f["/Header"]
        time = header.attrs.get("Time")
        time*= 0.9778e9 # to yr
        time/=1.e6 # to Myr
        xyz_p = np.array(f["/PartType0/Coordinates"])
        vel_p = np.array(f["/PartType0/Velocities"])
        mass_p = np.array(f["/PartType0/Masses"])
    rad2 = np.sum(xyz_p**2,axis=1)
    rad = np.sqrt(rad2)
#    zed = np.abs(xyz_p[:,2])
    selected = (rad<cutrad)# & (zed<cutzed)
    xyz_p = xyz_p[selected]
    vel_p = vel_p[selected]
    mass_p = mass_p[selected]
    rad = rad[selected]
    mtot = np.sum(mass_p)
    angular_momentum = np.sum(mass_p*np.sqrt(np.sum(np.cross(xyz_p,vel_p)**2,axis=1)))/mtot
    radial_momentum = np.sum(np.sum(xyz_p*vel_p,axis=1)/rad*mass_p)/mtot
#     geom_mean_rad = 1./(np.sum(1./(rad*mass_p))/mtot)
    minrad = np.min(rad)
    return time,angular_momentum,radial_momentum,mtot,minrad

if __name__ == '__main__':
    run_ids     = ["2020"]*5
    output_dirs = ["restest0m005_tiny","restest0m01_tiny","restest0m02_tiny","restest0m01_small","restest0m02_small"]
#     run_ids     = ["2020"]*2
#     output_dirs = ["restest0m01_small","restest0m02_small"]

    for run_id,output_dir in zip(run_ids,output_dirs):
        print("Dumping full evolution",run_id,output_dir)
        print("Reading and analysing angles")

        snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir)
        step = 1
        snapi = step
        snap_ids = range(snapi,snapf,step)
        with Pool(processes=128) as pool:
            output_data = pool.starmap(load_calc_summary,zip(it.repeat(run_id),it.repeat(output_dir),snap_ids))
        output_data=np.array(output_data)
        np.savetxt("data/dynamics{}_{}.dat".format(run_id,output_dir),output_data)