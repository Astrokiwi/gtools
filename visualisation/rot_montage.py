print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
#from joblib import Parallel, delayed

#from sys import path
#path.append("../src/")
import sph_frame

from sys import path
path.append("../")
import gizmo_tools

import matplotlib.pyplot as P

print("Running")

# output_dir = "q2redo"
# runs = ["q2edd05redo","q2edd10_aniso1redo","q2edd10_aniso3redo","q2edd10redo","q2edd20redo","q2redo"]
# runs = ["q2redo"]
# run_id = "2014"

runs = ["run_a2_e01"]
run_id = "2022"

snapx = 20
phis = np.linspace(0.,np.pi/2.,5) # representative
nsp = 5
# phis = np.linspace(0.,np.pi/2.,3) # representative
# nsp = 3
nphis = phis.size

gizmoDir = gizmo_tools.getGizmoDir(run_id)
movieDir = gizmo_tools.getMovieDir()

# rot = [0.,0.]
plot_thing = ['view']
L=32
width = .8
corners_face = [-width/2.,-width/2.]
corners_side = [0.,-width/2.]

view_cmap = "plasma"

for output_dir in runs:
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    infile = fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5"
    
    plotLabel, plotRanges, plotSliceTypes, plotCustomMass, plotData, logSliceTypes, extraBarTypes, plusMinusTypes, divergingTypes, customCmaps, customCmaps2 = sph_frame.pack_dicts()


    fig,sp = P.subplots(1,6,figsize=(12.,2.5), gridspec_kw = {'width_ratios':[1,16,16,16,16,16],'height_ratios':[16]})
#     fig,sp = P.subplots(1,4,figsize=(6.,2.2), gridspec_kw = {'width_ratios':[1,16,16,16],'height_ratios':[16]})
    cb_sp = sp[0]

    for icol in range(2,nsp+1):
        sp[icol].set_yticklabels([])
        sp[icol].get_shared_x_axes().join(sp[icol-1], sp[icol])
        sp[icol].get_shared_y_axes().join(sp[icol-1], sp[icol])

    for icol in range(1,nsp+1):
        sp[icol].set(adjustable='box-forced', aspect='equal')
#         sp[icol].set_xticks([-6,-4,-2,0,2,4,6])

    for icol,phi in enumerate(phis):
        rot = [0.,phi]
        time,data,x,y,z,rad2d,deep_face,deep_side,mask,n,n_ones = sph_frame.load_process_gadget_data(infile,rot,plot_thing,plotData,ringPlot=False,flatPlot=True)
        sph_frame.makesph_plot(fig,sp[icol+1],cb_sp,x,y,deep_face,0.,[data.brightness,data.opac,data.brightness,data.opac],data.m_p,data.h_p,L,mask,corners_face,width,r"$\log_{10} F$ (arbitrary units)",1.,7.,view_cmap,sph_frame.viewslice)
        sp[icol+1].set_title(r"$\phi=%2d^\circ$"%(phi*180./np.pi))
        
    sp[nsp].yaxis.set_label_position("right")
    sp[nsp].set_ylabel("t={0:.4f} Myr".format(time),size='x-large')

    cb_sp.yaxis.tick_left()
    cb_sp.yaxis.set_label_position("left")

    sp[1].set_xlabel("pc")
#     sp[1].set_ylabel("pc")
    
    fig.subplots_adjust(hspace=0., wspace=0.) 
    fig.tight_layout(pad=0.3,w_pad=0.0,h_pad=0.)
    P.savefig("../../figures/rotview_montage_"+run_id+output_dir+".png",dpi=150)
#     P.savefig("../../figures/rotview_montage_poster_"+run_id+output_dir+".png",dpi=150)
