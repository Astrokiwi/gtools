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

output_dir = "q2redo"
run_id = "2014"

snap_str = "1000"

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

ncuts = 3

gizmoDir = gizmo_tools.getGizmoDir()
movieDir = gizmo_tools.getMovieDir()

flattenedPlot = True
rot = [0.,0.]
plot_things = ['velr','TK_p','dustTemp']
load_things = ['vels','tdust','temp']
ranges = [[0.,2.5],[1.,4.],[1.,2.4]]
labels = [r'$\log v_r$ (km/s)',r'$\log T_g$ (K)',r'$\log T_d$ (K)']
# load_things = ['vels','temp']
L=256
width = 70.
corners_face = [-width/2.,-width/2.]
corners_side = [0.,-width/2.]

cmap = "plasma"

fullDir = gizmoDir+"/"+run_id+"/"+output_dir
infile = fullDir+"/snapshot_"+snap_str+".hdf5"


plotLabel, plotRanges, plotSliceTypes, plotCustomMass, plotData, logSliceTypes, extraBarTypes, plusMinusTypes, divergingTypes, customCmaps, customCmaps2 = sph_frame.pack_dicts()
time,data,x,y,z,rad2d,deep_face,deep_side,mask,n,n_ones = sph_frame.load_process_gadget_data(infile,rot,load_things,plotData,ringPlot=flattenedPlot,flatPlot=flattenedPlot)


fig,sp = P.subplots(ncuts,2,figsize=(4.,8.), gridspec_kw = {'width_ratios':[1,16],'height_ratios':[16]*ncuts})
# cb_sp = sp[0,0]

for irow in range(ncuts-1):
    for icol in [1]:
        sp[irow,icol].set_xticklabels([])
        sp[irow,icol].set_xticklabels([])
        sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow+1,icol])

for irow in range(ncuts):
    for icol in [1]:
        sp[irow,icol].set(adjustable='box-forced', aspect='equal')

for irow,plot_thing in enumerate(plot_things):
    v_p=data.__dict__[plot_thing]
    
    sph_frame.makesph_plot(fig,sp[irow,1],sp[irow,0],rad2d,z,deep_side,0.,v_p,data.m_p,data.h_p,L,mask,corners_side,width,labels[irow],ranges[irow][0],ranges[irow][1],cmap,sph_frame.weightslice)

for irow in range(ncuts):
    sp[irow,0].yaxis.tick_left()
    sp[irow,0].yaxis.set_label_position("left")

sp[ncuts-1,1].set_xlabel("pc")
sp[ncuts-1,1].set_ylabel("pc")

fig.subplots_adjust(hspace=0., wspace=1.) 
fig.tight_layout(pad=1.0,w_pad=0.0,h_pad=0.)
P.savefig("../../figures/outflow_montage_"+run_id+output_dir+".png",dpi=150)
