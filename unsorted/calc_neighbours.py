import numpy as np
import sys
import gizmo_tools
import h5py
# import scipy.spatial.distance as ds
sys.path.append("../src/")
from fort_calc_neighbours import fort_calc_neighbours

run_id = sys.argv[1]
output_dir = sys.argv[2]
snap_str = sys.argv[3]

gizmoDir = gizmo_tools.getGizmoDir(run_id)
with h5py.File(gizmoDir+"/"+run_id+"/"+output_dir+"/snapshot_"+snap_str+".hdf5","r") as f:
    xyz = np.array(f["/PartType0/Coordinates"])
    h = np.array(f["/PartType0/SmoothingLength"])
    print(np.min(h),np.max(h),np.min(xyz),np.max(xyz))
    neigh = fort_calc_neighbours.get_neigh(xyz,h)
#     neigh = np.zeros(h.size,dtype=int)
#     for ip,xyz_p in enumerate(xyz):
#         if ip%100==0: print(ip)
#         dists = np.sqrt(np.sum((xyz-xyz_p)**2,axis=1))
#         neigh[ip]=np.sum(dists<h)-1
    outp = np.array([xyz[:,0],xyz[:,1],xyz[:,2],neigh,h]).T
    np.savetxt("data/nneigh{}{}{}.dat".format(run_id,output_dir,snap_str),outp)
    