print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from matplotlib.mlab import PCA

import matplotlib.pyplot as P

import os
from sys import path
path.append("src/")
import tab_interp

from mpl_toolkits.mplot3d import Axes3D

print("Loading table (short form)")
chTab = tab_interp.CoolHeatTab(("coolheat_tab_marta/shrunk_table_labels_130617.dat"),("coolheat_tab_marta/shrunk_table_130617.dat"))
interpTabVec = np.vectorize(chTab.interpTab)

print("Running")

# luminosity = 6.2849e42 # erg/s

run_id = sys.argv[1]
output_dir = sys.argv[2]
snap_str = sys.argv[3]

f = h5py.File("/export/1/djw/gizmos/"+run_id+"/"+output_dir+"/snapshot_"+snap_str+".hdf5","r")

u_p = np.array(f["/PartType0/InternalEnergy"])
u_p*=1.e10 # to erg/g

rho_p = np.array(f["/PartType0/Density"])
rho_p*=6.77e-22 # to g/cm**3

# agn_heat_p = np.array(f["/PartType0/AGNHeat"])
# agn_heat_p*=1e10/3.08568e+16# to erg/s/g

flux_p = np.array(f["/PartType0/AGNIntensity"]) # energy per surface area per time
flux_p*=1.989e+53/(3.086e21)**2/(3.08568e+16)

depth_p = np.array(f["/PartType0/AGNColDens"]) # surface density units
depth_p*=(1.989e+43/3.086e+21**2) # to g/cm**2


molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16

depth_p/=(molecular_mass*proton_mass_cgs) # N in cm**(-2)
nH_p = rho_p/(molecular_mass*proton_mass_cgs)
TK_p = gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*u_p

tabStructs = interpTabVec(nH_p.astype(np.float64),TK_p.astype(np.float64),flux_p.astype(np.float64),depth_p.astype(np.float64))
heatRate = map(lambda y: y.dHeat, tabStructs)
heatRate = 10.**np.array(heatRate)
coolRate = map(lambda y: y.dCool, tabStructs)
coolRate = 10.**np.array(coolRate)

values = np.log10(np.array([nH_p,flux_p,depth_p,TK_p,heatRate])).T

print("PCA-ing")

pca_results = PCA(values)

N_skip = 100

# plot first two components
pca_x = pca_results.Y[:,0]
pca_xlabel = r"${:.3f}\log n_H+{:.3f}\log \Phi+{:.3f}\log N_H+{:.3f}\log T_g+{:.3f}\log H$".format(pca_results.Wt[0,0],pca_results.Wt[0,1],pca_results.Wt[0,2],pca_results.Wt[0,3],pca_results.Wt[0,4])

pca_y = pca_results.Y[:,1]
pca_ylabel = r"${:.3f}\log n_H+{:.3f}\log \Phi+{:.3f}\log N_H+{:.3f}\log T_g+{:.3f}\log H$".format(pca_results.Wt[1,0],pca_results.Wt[1,1],pca_results.Wt[1,2],pca_results.Wt[1,3],pca_results.Wt[1,4])

pca_z = pca_results.Y[:,2]
pca_zlabel = r"${:.3f}\log n_H+{:.3f}\log \Phi+{:.3f}\log N_H+{:.3f}\log T_g+{:.3f}\log H$".format(pca_results.Wt[2,0],pca_results.Wt[2,1],pca_results.Wt[2,2],pca_results.Wt[2,3],pca_results.Wt[2,4])

# H,xedges,yedges = np.histogram2d(pca_x,pca_y,bins=(150,150))
# H = np.log10(H.T)

nbins = 25

H, edges = np.histogramdd(pca_results.Y[:,0:3], bins = (nbins,nbins,nbins))
nonzeroH = H!=0
nonzeroHpoints = np.where(H!=0)
nonzeroHdense = H[nonzeroH]

print("Plotting")

P.close()
fig = P.figure()
ax = fig.add_subplot(111, projection='3d')

#ax.scatter(nonzeroHpoints[0],nonzeroHpoints[1],nonzeroHpoints[2],s=((np.log10(nonzeroHdense)+1)**2)*16.,c='red',marker='x')
scat=ax.scatter(nonzeroHpoints[0],nonzeroHpoints[1],nonzeroHpoints[2],c=np.log10(nonzeroHdense),marker='x',cmap='plasma')
P.colorbar(scat,label=r"log count")

#P.pcolormesh(xedges,yedges,H,cmap='plasma',vmin=-1.,vmax=np.max(H))
# ax.scatter(pca_x,pca_y,pca_z)
# ax.set_xlabel(pca_xlabel)
# ax.set_ylabel(pca_ylabel)
# ax.set_zlabel(pca_zlabel)
# 
ax.set_xlabel("First component")
ax.set_ylabel("Second component")
ax.set_zlabel("Third component")

print("First component:"+pca_xlabel)
print("Second component:"+pca_ylabel)
print("Third component:"+pca_zlabel)

fig.suptitle("1st:"+pca_xlabel+"\n2nd:"+pca_ylabel+"\n3rd:"+pca_zlabel)

# # P.xlim(xedges[0],xedges[-1])
# # P.ylim(yedges[0],yedges[-1])
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1.))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(1.))
ax.zaxis.set_minor_locator(ticker.MultipleLocator(1.))
# ax.tick_params(which='both',direction="out")
#P.colorbar(label=r"$\log$ count")
for irot,angle in enumerate(np.linspace(0.,360.,60,endpoint=False)):
    fname = "pics/pca_heat_"+run_id+output_dir+snap_str+"{:03d}.png".format(irot)
    print("Rendering for plot, angle={} deg".format(angle))
    ax.view_init(30, angle)
    P.savefig(fname,dpi=200)

cmd = "ffmpeg -y -r 24 -i pics/pca_heat_"+run_id+output_dir+snap_str+"%03d.png -c:v mpeg4 -q:v 1 /export/1/djw/movies/PCArot"+run_id+output_dir+snap_str+".mp4"
print(cmd)
os.system(cmd)