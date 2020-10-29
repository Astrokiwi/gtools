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
path.append("../visualisation/")
import gizmo_tools

import matplotlib.pyplot as P

print("Running")

# output_dir = "q2redo"
# run_id = "2014"
# 
# snap_str = "1000"

# output_dirs = ["newtable00001","newtable00005","newtable","newtable0005","newtable001","newtable005","newtable01"]
# run_id = "2022"
# 
# snap_str = "039"
# 
# output_dirs = ["newtable00001","newtable00005","newtable","newtable0005","newtable001","newtable005","newtable01"]
# run_id = "2022"
# 
# snap_str = "039"

output_dirs = ["a2_e01_0m00001","a2_e01","a2_e01_0m001","a2_e01_0m01","a2_e01_0m1"]
run_id = "3001"

snap_str = "020"

# cuts = [    ['temp',10.**0.,10.**2.5],
#             ['temp',10.**2.5,10.**4.],
#             ['tdust',1.,30.],
#             ['tdust',30.,150.]
#             ]
# cuts = [    ['temp',10.**0.,10.**2.5],
#             ['temp',10.**2.5,10.**4.],
#             ['temp',10.**0.,10.**2.5],
#             ['temp',10.**2.5,10.**4.]
#             ]

# labels = [  r'$\log T_g\leq2.5$ (K)',
#             r'$\log T_g\geq2.5$ (K)',
#             r'$T_d\leq30.$ K',
#             r'$T_d\geq30.$ K'
#             ]


nruns = len(output_dirs)

gizmoDir = gizmo_tools.getGizmoDir(run_id)
movieDir = gizmo_tools.getMovieDir()

flattenedPlot = True
rot = [0.,0.]
load_things = ['nH']
# load_things = ['vels','temp']
L=256
width = 0.3
corners_face = [-width/2.,-width/2.]
corners_side = [0.,-width/2.]

cmap = "inferno"
# cmap = "bwr"


plotLabel, plotRanges, plotSliceTypes, plotCustomMass, plotData, logSliceTypes, extraBarTypes, plusMinusTypes, divergingTypes, customCmaps, customCmaps2 = sph_frame.pack_dicts()

nx = 2
ny = nruns//2
if ny*2<nruns: ny+=1

fig,sp = P.subplots(ny,1+nx,figsize=(8.,9.5), gridspec_kw = {'width_ratios':[16]*nx+[1],'height_ratios':[16]*ny})
# cb_sp = sp[0,0]

for irow in range(ny-1):
    for icol in [0,1]:
        sp[irow,icol].set_xticklabels([])
        sp[irow,icol].set_xticklabels([])
    for icol in [1]:
        sp[irow,icol].set_yticklabels([])
    for icol in [0]:
        sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow+1,icol])
        sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow,icol+1])

for irow in range(ny):
    for icol in [0,1]:
        sp[irow,icol].set(adjustable='box-forced', aspect='equal')

for irow,output_dir in enumerate(output_dirs):
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    infile = fullDir+"/snapshot_"+snap_str+".hdf5"
    time,data,x,y,z,rad2d,deep_face,deep_side,mask,n,n_ones = sph_frame.load_process_gadget_data(infile,rot,load_things,plotData,ringPlot=flattenedPlot,flatPlot=flattenedPlot)
    v_p=data.nH_p
    
#     ix = irow%2+1
    ix = irow%2
    iy = irow//2
    sph_frame.makesph_plot(fig,sp[iy,ix],sp[iy,2],rad2d,z,deep_side,0.,v_p,data.m_p,data.h_p,L,mask,corners_side,width,r"$\log n_H$ (cm$^{-3}$)",1.,9.,cmap,sph_frame.weightslice)

if nx*ny>nruns:
    sp[ny-1,nx].remove()
    sp[ny-1,nx-1].remove()

# for irow in range(ncuts):
#     sp[irow,0].yaxis.tick_left()
#     sp[irow,0].yaxis.set_label_position("left")

sp[ny-1,0].set_xlabel("pc")
sp[ny-1,0].set_ylabel("pc")

# fig.subplots_adjust(hspace=0., wspace=0.) 
fig.subplots_adjust(left=0.1,right=.9,bottom=0.05,top=.99,hspace=0.03,wspace=0.01) 
# fig.tight_layout(pad=0.0,w_pad=0.0,h_pad=0.)
# P.savefig("../../figures/res_montage_"+run_id+output_dir+".png",dpi=150)
P.savefig("../../figures/new_res_montage_"+run_id+output_dir+".png",dpi=150)
