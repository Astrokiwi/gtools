import sys
import gizmo_tools
import numpy as np
import h5py

molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16


if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    snapi = 100
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir)
    step = 1

    oldSnapshot = None
    
    max_dTK_p = 0
    max_dTK_id = None
    
    for isnap in range(snapi,snapf+1,step):
        snap_str = "%03d"%isnap
        gizmoDir = gizmo_tools.getGizmoDir()
        fullDir = gizmoDir+"/"+run_id+"/"+output_dir
        f = h5py.File(fullDir+"/snapshot_"+snap_str+".hdf5","r")
        thisSnapshot = dict()
        TK_in = gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*1.e10*np.array(f["/PartType0/InternalEnergy"])
        N = TK_in.size
        thisSnapshot["id_p"] = np.array(f["/PartType0/ParticleIDs"]).astype(int)
        thisSnapshot["TK_p"] = np.zeros(N)
        thisSnapshot["TK_p"][thisSnapshot["id_p"]-1] = TK_in
#         thisSnapshot["TK_p"] = TK_in

        if oldSnapshot:
            dTK_p = np.maximum(np.abs((oldSnapshot["TK_p"]-thisSnapshot["TK_p"])/oldSnapshot["TK_p"]),np.abs((oldSnapshot["TK_p"]-thisSnapshot["TK_p"])/thisSnapshot["TK_p"]))
            this_max_dTK_index = np.argmax(dTK_p)
            this_max_dTK_p = dTK_p[this_max_dTK_index]
            if this_max_dTK_p>max_dTK_p:
                max_dTK_p=this_max_dTK_p
                max_dTK_id=this_max_dTK_index+1
                print(snap_str," New max:",max_dTK_p,max_dTK_id,oldSnapshot["TK_p"][this_max_dTK_index],thisSnapshot["TK_p"][this_max_dTK_index])
            
        oldSnapshot = thisSnapshot
        
        
