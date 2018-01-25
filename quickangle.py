import sys
import gizmo_tools
import numpy as np
import h5py

bins = np.arange(0.,91.,.25)


if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    snap_str = sys.argv[3]

    gizmoDir = gizmo_tools.getGizmoDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    f = h5py.File(fullDir+"/snapshot_"+snap_str+".hdf5","r")

    header = f["/Header"]
    time = header.attrs.get("Time")
    time*= 0.9778e9 # to yr
    time/=1.e6 # to Myr

    xyz_p = np.array(f["/PartType0/Coordinates"])
    mass_p = np.array(f["/PartType0/Masses"])

    mass_p*=1.e10 # to Msun
    
    rad2d_p = np.sqrt(xyz_p[:,0]**2+xyz_p[:,1]**2)
    z_p = xyz_p[:,2]
    
    theta_p = np.abs(np.arctan2(z_p,rad2d_p)*180./np.pi)
    
    print(np.median(theta_p),np.mean(theta_p))
        
    theta_histogram, theta_edges = np.histogram(theta_p,weights=mass_p,bins=bins)
    
    output_data = np.array([theta_edges[:-1],theta_histogram]).T
    
    np.savetxt("data/quickangle"+run_id+output_dir+snap_str+".dat",output_data)
