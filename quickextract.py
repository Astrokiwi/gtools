import sys
import numpy as np
import h5py
import gizmo_tools

run_id = sys.argv[1]
output_dir = sys.argv[2]
snap_id = sys.argv[3]


gizmoDir = gizmo_tools.getGizmoDir()
fullDir = gizmoDir+"/"+run_id+"/"+output_dir
with h5py.File(fullDir+"/snapshot_"+snap_id+".hdf5","r") as f:
    h_p = np.array(f["/PartType0/SmoothingLength"])
    h_sort = np.sort(h_p)
    print(h_p[-20:])