import sys
import numpy as np
import h5py
import gizmo_tools
# import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as P


run_id = sys.argv[1]
output_dir = sys.argv[2]
snap_id = sys.argv[3]


gizmoDir = gizmo_tools.getGizmoDir(run_id)
fullDir = gizmoDir+"/"+run_id+"/"+output_dir
with h5py.File(fullDir+"/snapshot_"+snap_id+".hdf5","r") as f:
#     v = np.array(f["/PartType0/TimeStep"])
#     v = np.array(f["/PartType0/RadiativeAcceleration"])
#     v = np.array(f["/PartType0/AGNHeat"])
    v2 = np.array(f["/PartType0/TimeStep"])
    v = np.array(f["/PartType0/InternalEnergy"])
    v = v[v2>1.e-9]
    
    vsort = np.sort(v)
    print(vsort[-20:])
    print(vsort[:20])
    print(np.min(v),np.mean(v),np.median(v),np.max(v))
    
#     np.savetxt("data/quickdump.dat",[v,v2])
#     heatTime = v/v2
#     heatTimeSorted = np.sort(heatTime)
#     print(heatTimeSorted[140179],heatTimeSorted[-140179])
#     print(np.sum(heatTimeSorted<0.),np.sum(heatTimeSorted>0.),np.sum(np.abs(heatTimeSorted)<1.e-9))

    
#     heat = np.array(f["/PartType0/AGNHeat"])
#     heat_sort = np.sort(heat)
#     print(heat_sort[-20:],np.max(heat),np.min(heat),np.sum(heat!=0.))