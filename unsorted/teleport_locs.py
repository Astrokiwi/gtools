import sys
import gizmo_tools
import numpy as np
import h5py

    
if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    
    snapi = 0
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir)

    gizmoDir = gizmo_tools.getGizmoDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir

    output = np.empty((0,2))
    print(output.shape)
    
#    snapf = 100

    for i in range(snapf):
        f = h5py.File(fullDir+"/snapshot_"+"%03d" % i+".hdf5","r")

        header = f["/Header"]
        time = header.attrs.get("Time")
        time*= 0.9778e9 # to yr
        time/=1.e6 # to Myr

        xyz_p = np.array(f["/PartType0/Coordinates"])
        vel_p = np.array(f["/PartType0/Velocities"])
        
        id_p = np.array(f["/PartType0/ParticleIDs"])
        
        # check I'm sorting this right
        xyz_p[id_p-1,:] = xyz_p.copy()
        vel_p[id_p-1,:] = vel_p.copy()
        
        if i!=0:
            dt = time-old_time
            d_xyz = dt * old_vel_p * 1.02269032e-3 # unit conversion from 1 Myr km/s to kpc

            predict_xyz = old_xyz_p + d_xyz
            error2_xyz = np.sum((predict_xyz - xyz_p)**2,axis=1)
            d2_xyz = np.sum(d_xyz**2,axis=1)
            max_error_xyz = error2_xyz/d2_xyz
#             print(old_xyz_p[0,:])
#             print(xyz_p[0,:])
#             print((xyz_p-old_xyz_p)[0,:]/d_xyz[0,:])
#             print((xyz_p-old_xyz_p)[0,:]/dt)
#             print(d_xyz[0,:])
#             print(old_vel_p[0,:])
#             print(vel_p[0,:])
#             print(d_xyz)
#             print(predict_xyz)
#             print(xyz_p)
#             print(error2_xyz)
#             print(d2_xyz)
#             print(max_error_xyz)
#             sys.exit()
            teleported = (max_error_xyz>10000.)
            n_error = np.sum(teleported)
            
            output_to_append = np.full((n_error,2),time)
            output_to_append[:,1] = np.sqrt(np.sum(old_xyz_p[teleported,:]**2,axis=1))
            
            output = np.vstack([output,output_to_append])
            print(i,time,n_error,output.shape)
        
        old_xyz_p = xyz_p
        old_vel_p = vel_p
        old_time = time
    np.savetxt("data/sf_locs"+run_id+output_dir+".dat",output)




