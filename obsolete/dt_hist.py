import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import h5py
import sys
import matplotlib.pyplot as P

nruns = (len(sys.argv)-1)/3

fig = P.figure()

inbins = np.linspace(1.6,3.,12)

alldt = []
alllabels = []

for irun in range(nruns):

    run_id = sys.argv[1+irun*3]
    output_dir = sys.argv[2+irun*3]
    snap_str = sys.argv[3+irun*3]

    f = h5py.File("/export/1/djw/gizmos/"+run_id+"/"+output_dir+"/snapshot_"+snap_str+".hdf5","r")
    
    dt_p = np.array(f["/PartType0/TimeStep"])
    dt_p*=0.9778e9 # to yr
    dt_p = np.log10(dt_p)
    
    alldt.append(dt_p)
    alllabels.append(output_dir)

n, bins, patches = P.hist(alldt,bins=inbins,label=alllabels,alpha=.5,histtype='step')

P.legend()
P.savefig("pics/manydthist.png",dpi=300)
P.close()