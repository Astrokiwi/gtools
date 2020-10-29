print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys
import matplotlib as mpl
mpl.use('Agg')

import pylab as P

print("Running")


run_id = sys.argv[1]
output_dir = sys.argv[2]
snap_str = sys.argv[3]

f = h5py.File("/export/1/djw/gizmos/"+run_id+"/"+output_dir+"/snapshot_"+snap_str+".hdf5","r")

xyz = np.array(f["/PartType0/Coordinates"])
rad2d_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
phi = np.abs(np.arctan2(xyz[:,2],rad2d_p))/(2.*np.pi)*360

m_p = np.array(f["/PartType0/Masses"])
#m_p*=1.989e+43 # 10^10 solar masses to g
m_p*=1.e10 # 10^10 solar masses to g

phibins = np.linspace(0.,90.,19)

P.figure()
n, bins, patches = P.hist(phi, bins=phibins,weights=m_p,color='white',log=True)
P.xlabel(r"$\phi$ ($^\degree$)")
P.ylabel(r"$M$ ($M_\odot$)")
P.savefig("pics/m_angle_"+run_id+"_"+output_dir+"_"+snap_str+".png",dpi=200)
P.close()
