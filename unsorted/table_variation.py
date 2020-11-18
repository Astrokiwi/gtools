print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.colors as colors
import matplotlib.ticker as ticker

import matplotlib.pyplot as P

from scipy import stats

from sys import path
path.append("../src/")
import tab_interp

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

values = [nH_p,flux_p,depth_p,TK_p]
#values = [nH_p,flux_p,depth_p,TK_p]
titles = ["nH","flux","col","Tg"]
labels = [r"$\log n_H$ (cm$^{-3}$)",r"$\log \Phi$ (erg/cm$^{-2}$/s)",r"$\log N_H$ (cm$^{-2}$)",r"$\log T_g$ (K)"]

includedVals = [[0,3],[1,2]]

heatCoolNames = ["heat","cool"]

# if ( not depth_exists ):
#     nv = len(values)-1
# else:
#     nv = len(values)
nv = len(values)


for ivals in includedVals:
    iv = ivals[0]
    jv = ivals[1]
    for stat in 'mean','median','max','min':
        P.figure()
        vx = np.log10(values[iv])
        vy = np.log10(values[jv])
        notnans = ( (np.isfinite(vx)) & (np.isfinite(vy)) & (np.isfinite(heatRate)) )
        vx = vx[notnans]
        vy = vy[notnans]
        vz = heatRate[notnans]
        H,xedges,yedges,nbins = stats.binned_statistic_2d(vx,vy,vz,statistic=stat,bins=(50,50))
        print(H)
        #H_norm,xedges,yedges = np.histogram2d(vy,vx,bins=(150,150))
        #H = H/H_norm
        H = np.log10(H).T
        #H = H.T
        P.pcolormesh(xedges,yedges,H,cmap='plasma',vmin=-23.,vmax=-10.) #,norm=colors.LogNorm()
    
        P.xlabel(labels[iv])
        P.ylabel(labels[jv])
        P.xlim(xedges[0],xedges[-1])
        P.ylim(yedges[0],yedges[-1])
        ax = P.gca()
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(1.))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(1.))
        P.tick_params(which='both',direction="out")
        #P.colorbar(label=r"$\log$ count")
        P.colorbar(label=r"$\log H$")
        fname = "../figures/"+stat+run_id+output_dir+titles[iv]+titles[jv]+snap_str+".png"
        P.savefig(fname,dpi=150)
        P.close()



