import numpy as np
import io
import sys

import gizmo_tools


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# import matplotlib as mpl
# import matplotlib.pyplot as P

t_in_yr = 0.9778e9

def dump_performance(run_id,output_dir):
    print(run_id,output_dir)
    gizmoDir = gizmo_tools.getGizmoDir(run_id)
    movieDir = gizmo_tools.getMovieDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir


    filename = fullDir+"/info.txt"

    s = open(filename).read().replace(',',' ')

    data = np.loadtxt(io.StringIO(s),dtype=str)

    sync_point = data[:,1].astype(int)

    time = data[:,3].astype(float) * t_in_yr

    n_active = data[:,6].astype(float)

    dt = data[:,8].astype(float) * t_in_yr

    n_smooth = moving_average(n_active,n=4096)
    dt_smooth = moving_average(dt,n=4096)

    dn = n_active.size-n_smooth.size
    dn_before = dn//2
    dn_after = dn-dn_before

    n_smooth = np.lib.pad(n_smooth,(dn_before,dn_after),'constant',constant_values=(0,0))
    dt_smooth = np.lib.pad(dt_smooth,(dn_before,dn_after),'constant',constant_values=(0,0))

    out_data = np.array([sync_point,time,n_active,dt,n_smooth,dt_smooth]).T

    np.savetxt("data/info_plotable_"+run_id+output_dir+".dat",out_data)

if __name__== "__main__":
    run_id = sys.argv[1]
    output_dir = sys.argv[2]

    dump_performance(run_id,output_dir)
