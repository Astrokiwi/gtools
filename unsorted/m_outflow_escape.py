import sys
import gizmo_tools
import numpy as np
import h5py

G_kms_pc_msun = 0.0043022682

def v_esc(r,m_bh,m_hern,a_hern):
    return np.sqrt(2*G_kms_pc_msun*(m_bh/r + m_hern/(a_hern+r)))

def m_outflow(run_id,output_dir,snap_str):
    gizmoDir = gizmo_tools.getGizmoDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    f = h5py.File(fullDir+"/snapshot_"+snap_str+".hdf5","r")

    header = f["/Header"]
    time = header.attrs.get("Time")
    time*= 0.9778e9 # to yr
    time/=1.e6 # to Myr

    xyz_p = np.array(f["/PartType0/Coordinates"])
    vel_p = np.array(f["/PartType0/Velocities"])
    mass_p = np.array(f["/PartType0/Masses"])

    mass_p*=1.e10 # to Msun
    xyz_p*=1.e3 # to pc
    
    rad3d_p = np.sqrt(np.sum(xyz_p**2,axis=1))
    vmag = np.sqrt(np.sum(vel_p**2,axis=1))
    
    outflowing = (v_esc(rad3d_p,1.e6,1.e9,250.)<vmag)
    
    m_outflowing = np.sum(mass_p[outflowing])
    
    return time,m_outflowing


if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]

    snapi = 0
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir)
    snap_strs = ["%03d" % i for i in range(snapf+1)]
    with Pool(processes=80) as pool:
        output_data = pool.starmap(m_outflow,zip(it.repeat(run_id),it.repeat(output_dir),snap_strs))
    np.savetxt("data/m_outflow"+run_id+output_dir+".dat",output_data)
