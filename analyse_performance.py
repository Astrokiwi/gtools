import numpy as np
import StringIO
import sys

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# import matplotlib as mpl
# import matplotlib.pyplot as P

t_in_yr = 0.9778e9

run_id = sys.argv[1]
output_dir = sys.argv[2]

filename = "/export/1/djw/gizmos/"+run_id+"/"+output_dir+"/info.txt"

s = open(filename).read().replace(',',' ')

data = np.loadtxt(StringIO.StringIO(s),dtype=str)

sync_point = data[:,1].astype(int)

time = data[:,3].astype(float) * t_in_yr

n_active = data[:,6].astype(float)

n_smooth = moving_average(n_active,n=5000)

dn = n_active.size-n_smooth.size
dn_before = dn/2
dn_after = dn-dn_before

n_smooth = np.lib.pad(n_smooth,(dn_before,dn_after),'constant',constant_values=(0,0))

dt = data[:,8].astype(float) * t_in_yr

out_data = np.array([sync_point,time,n_active,dt,n_smooth]).T

np.savetxt("data/info_plotable_"+run_id+output_dir+".dat",out_data)