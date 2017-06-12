print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.colors as colors

import pylab as P

from sys import path
path.append("src/")
import tab_interp

print("Loading table (short form)")
chTab = tab_interp.CoolHeatTab(("coolheat_tab_marta/shrunk_table_labels_080517.dat"),("coolheat_tab_marta/shrunk_table_080517.dat"))
interpTabVec = np.vectorize(chTab.interpTab)

print("Running")

# luminosity = 6.2849e42 # erg/s

run_id = sys.argv[1]
output_dir = sys.argv[2]
snap_str = sys.argv[3]

f = h5py.File("/export/1/djw/gizmos/"+run_id+"/"+output_dir+"/snapshot_"+snap_str+".hdf5","r")

xyz = np.array(f["/PartType0/Coordinates"])
rad_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2+xyz[:,2]**2)

rad2d_p = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
rad2d_p*=1. # already in kpc

z_p = np.abs(xyz[:,2])
z_p*=1. # already in kpc

u_p = np.array(f["/PartType0/InternalEnergy"])
u_p*=1.e10 # to erg/g

rho_p = np.array(f["/PartType0/Density"])
rho_p*=6.77e-22 # to g/cm**3

agn_heat_p = np.array(f["/PartType0/AGNHeat"])
agn_heat_p*=1e10/3.08568e+16# to erg/s/g

dt_p = np.array(f["/PartType0/TimeStep"])


N_part = rho_p.size

# if "/PartType0/AGNColDens" in f:
#     depth_p = np.array(f["/PartType0/AGNColDens"]) # surface density units
#     depth_exists = True
# else:
#     depth_p = np.zeros(N_part)
#     depth_exists = False

flux_p = np.array(f["/PartType0/AGNIntensity"]) # energy per surface area per time
flux_p*=1.989e+53/(3.086e21)**2/(3.08568e+16)

depth_p = np.array(f["/PartType0/AGNColDens"]) # surface density units
depth_p*=(1.989e+43/3.086e+21**2) # to g/cm**2

radaccel_p = np.array(f["/PartType0/RadiativeAcceleration"])
radaccel_p*=3.0857e21/3.08568e+16**2 # to cm/s/s

radrad_p = (xyz[:,0]*radaccel_p[:,0]+xyz[:,1]*radaccel_p[:,1]+xyz[:,2]*radaccel_p[:,2])/rad_p

molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16

depth_p/=(molecular_mass*proton_mass_cgs) # N in cm**(-2)
# flux_p = luminosity/(4.*np.pi*rad_p**2)
# flux_p/=(3.086e21)**2 # erg/s/kpc**2 to erg/s/cm**2

nH_p = rho_p/(molecular_mass*proton_mass_cgs)
TK_p = gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*u_p

tabStructs = interpTabVec(nH_p.astype(np.float64),TK_p.astype(np.float64),flux_p.astype(np.float64),depth_p.astype(np.float64))
dustTemp = map(lambda y: y.dustT, tabStructs)
dustTemp = np.array(dustTemp)

#nH_p = 10.**((4.-np.log10(TK_p))**2)

# values = [nH_p,flux_p,TK_p,depth_p]
# titles = ["nH","flux","TK","col"]
# labels = [r"$\log n_H$ (cm$^{-3}$)",r"$\log \Phi$ (erg/cm$^{-2}$/s)",r"$\log T$ (K)",r"$\log N_H$ (cm$^{-2}$)"]

values = [nH_p,flux_p,depth_p,TK_p,rad2d_p,z_p,dustTemp,radrad_p]
#values = [nH_p,flux_p,depth_p,TK_p]
titles = ["nH","flux","col","Tg","r","z","Td","arad"]
labels = [r"$\log n_H$ (cm$^{-3}$)",r"$\log \Phi$ (erg/cm$^{-2}$/s)",r"$\log N_H$ (cm$^{-2}$)",r"$\log T_g$ (K)",r"$r$ (kpc)",r"$|z|$ (kpc)",r"$T_d$ (K)",r"$a_{rad}$ (cm/s/s)"]
dolog = [True,True,True,True,False,False,False,True]
ranges = [None,None,None,None,[0,.02],[0.,.02],None,None]
#ranges = [None,None,None,None,None,None,None,None]

includedVals = [0,1,2,3,4,5,6,7]

# if ( not depth_exists ):
#     nv = len(values)-1
# else:
#     nv = len(values)
nv = len(values)

#bigslice = (dustTemp>235.) & (dustTemp<245.) & (depth_p<1.e15)

#bigslice = (TK_p>1.e6)

#bigslice = (agn_heat_p<-1.e6)

#bigslice = (dt_p<1.e-10)

bigslice = (dt_p>0.)
#bigslice = (TK_p>10.**3.5) & (TK_p<10.**4.) & (depth_p>10.**22.) & (depth_p<10.**23.2)

# for i,iv in enumerate(includedVals):
#     for jv in includedVals[i+1:]:
for iv in range(7):
    for jv in [7]:
        fname = "../figures/"+run_id+output_dir+titles[iv]+titles[jv]+snap_str+".png"
        P.figure()
        P.xlabel(labels[iv])
        P.ylabel(labels[jv])
        if ( dolog[iv] ):
            vx = np.log10(values[iv])
        else:
            vx = values[iv]
        if ( dolog[jv] ):
            vy = np.log10(values[jv])
        else:
            vy = values[jv]
        notnans = ((np.isfinite(vx)) & (np.isfinite(vy))) & bigslice
        if ( not (ranges[iv] is None) ):
            notnans = notnans & (vx>=ranges[iv][0]) & (vx<=ranges[iv][1])
        if ( not (ranges[jv] is None) ):
            notnans = notnans & (vy>=ranges[jv][0]) & (vy<=ranges[jv][1])
        vx = vx[notnans]
        vy = vy[notnans]
        H,xedges,yedges = np.histogram2d(vx,vy,bins=(150,150))
        H = np.log10(H).T
        #H = H.T
        P.pcolormesh(xedges,yedges,H,cmap='plasma',vmin=np.min(H[np.isfinite(H)]),vmax=np.max(H[np.isfinite(H)])) #,norm=colors.LogNorm()
        P.xlim(xedges[0],xedges[-1])
        P.ylim(yedges[0],yedges[-1])
        #P.colorbar(label=r"$\log$ count")
        P.colorbar(label=r"count")
        P.savefig(fname,dpi=150)
        P.close()



