print("Importing")

# Import libraries to do our magic
import numpy as np
import h5py
import sys
import matplotlib as mpl
mpl.use('Agg')

import pylab as P

import os
#from joblib import Parallel, delayed

#from sys import path
#path.append("../src/")
from sph_plotter import sph_plotter

print("Running")


run_id = sys.argv[1]
output_dir = sys.argv[2]

fullDir = "/export/1/djw/gizmos/"+run_id+"/"+output_dir


cmaps = ["viridis","inferno","plasma","magma","jet","Greys","cool","winter","hot","gray","terrain","ocean","rainbow","gnuplot","Accent","seismic"]

lfiggrid = int(np.ceil(np.sqrt(len(cmaps))))

#for snapx in range(2,20):
for snapx in range(10,11):
    print("Loading",snapx)
    f = h5py.File(fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5","r")


    xyz = np.array(f["/PartType0/Coordinates"])
    m_p = np.array(f["/PartType0/Masses"])
    h_p = np.array(f["/PartType0/SmoothingLength"])
    u_p = np.array(f["/PartType0/InternalEnergy"])

    u_p*=1.e10 # to erg/g
    m_p*=1e+10 # 10^10 solar masses to solar masses
    molecular_mass = 4./(1.+3.*.76)
    proton_mass_cgs = 1.6726e-24
    gamma_minus_one = 5./3.-1.
    boltzmann_cgs = 1.38066e-16

    TK_p = (gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)*u_p)


    n = h_p.size
    #L = 256
    L = 512
    x = xyz[:,0]*1.e3
    y = xyz[:,1]*1.e3
    z = xyz[:,2]*1.e3
    h_p*=1.e3

    mask = np.full(n,True,dtype=bool)
    n_ones = np.ones(n)


    fig, sp = P.subplots(lfiggrid,lfiggrid)
    #fig.suptitle(r"$T="+("%.4f" % time)+"$ Myr")
    fw_inches = 12.
    fig.set_figwidth(fw_inches)
    fig.set_figheight(fw_inches)

    corner = [-10.,-10.]
    width = 20.

    print("SPH-ing")
    map = sph_plotter.sph_weight(x,y,m_p,h_p,TK_p,L,corner,width,mask,n)
    map = np.log10(map)

    print(np.max(map),np.min(map))

    xedges = np.arange(corner[0],corner[0]+width,width/L)
    yedges = np.arange(corner[1],corner[1]+width,width/L)
    zlow = 2.5
    zhigh = 7.5

    for icmap in range(len(cmaps)):
        print("Plotting %d"%icmap)
        this_sp = sp[icmap//lfiggrid,icmap%lfiggrid]
        cmap_label = cmaps[icmap]
        this_cmap = P.get_cmap(cmap_label)
        this_cmap.set_bad('black',1.)
        #this_cmap.set_under('black',1.)
        mesh = this_sp.pcolormesh(xedges,yedges,map.T,cmap=this_cmap,vmin=zlow,vmax=zhigh) 
        #mesh = sp.pcolormesh(xedges,yedges,map.T,cmap=this_cmap) 
        this_sp.axis('equal') 
        this_sp.set_title(cmap_label)
    
    
    pos = sp[0,0].get_position()
    lpixx = (L/(pos.x1-pos.x0))
    my_dpi = int(np.floor(lpixx/fw_inches))

    print("Rendering & saving")
    outfile = "../pics/sph_colourtest_big%03d.png"%snapx
    #fig.subplots_adjust(left=0.07,hspace=.07,bottom=.1,top=.9)
    fig.tight_layout()
    #fig.savefig(outfile,dpi=200)
    fig.savefig(outfile,dpi=my_dpi)
