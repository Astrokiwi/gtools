print("Importing")

# Import libraries to do our magic

from ..tools import sph_frame
from ..tools import gizmo_tools

import matplotlib.pyplot as plt

print("Running")

# output_dir = "q2redo"
# run_id = "2014"
# 
# snap_str = "1000"

# output_dir = "run_a2_e01"
# run_id = "2022"
output_dir = "a2_e01"
run_id = "3001"

snap_str = "020"


# for paper
# cuts = [    ['temp',10.**0.,10.**3.],
#             ['temp',10.**3.,10.**5.],
#             ['tdust',10.,10.**2.],
#             ['tdust',10.**2.,10.**4]
#             ]
# labels = [  r'$\log T_g\leq3$ (K)',
#             r'$\log T_g\geq3$ (K)',
#             r'$\log T_d\leq2$ (K)',
#             r'$\log T_d\geq2$ (K)'
#             ]

# for TORUS2018
cuts = [    ['tdust',10.,10.**2.],
            ['tdust',10.**2.,10.**4]
            ]
labels = [  r'$\log T_d\leq2$ (K)',
            r'$\log T_d\geq2$ (K)'
            ]

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

ncuts = len(cuts)

gizmoDir = gizmo_tools.getGizmoDir(run_id)
movieDir = gizmo_tools.getMovieDir()

flattenedPlot = True
rot = [0.,0.]
plot_thing = ['vels']
load_things = ['vels','tdust','temp']
# load_things = ['vels','temp']
L=256
width = .8
# width = .6
corners_face = [-width/2.,-width/2.]
corners_side = [0.,-width/2.]

cmap = "plasma"

fullDir = gizmoDir+"/"+run_id+"/"+output_dir
infile = fullDir+"/snapshot_"+snap_str+".hdf5"


plotLabel, plotRanges, plotSliceTypes, plotCustomMass, plotData, logSliceTypes, extraBarTypes, plusMinusTypes, divergingTypes, customCmaps, customCmaps2 = sph_frame.pack_dicts()
time,data,x,y,z,rad2d,deep_face,deep_side,mask_default,n,n_ones = sph_frame.load_process_gadget_data(infile,rot,load_things,plotData,ringPlot=flattenedPlot,flatPlot=flattenedPlot)

print("time=",time)

fig,sp = plt.subplots(ncuts, 3, figsize=(6., 2.5 * ncuts), gridspec_kw = {'width_ratios':[1, 16, 16], 'height_ratios': [16] * ncuts})
cb_sp = sp[0,0]

for irow in range(ncuts-1):
    for icol in [1,2]:
        sp[irow,icol].set_xticklabels([])
        sp[irow,icol].set_xticklabels([])
        sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow+1,icol])

for irow in range(ncuts):
    for icol in [2]:
        sp[irow,icol].set_yticklabels([])
    for icol in [1,2]:
        sp[irow,icol].set(adjustable='box-forced', aspect='equal')

for irow,cut in enumerate(cuts):
    v=data.__dict__[plotData[cut[0]]]
    mask = (v>cut[1]) & (v<cut[2])
    
    sph_frame.makesph_plot(fig,sp[irow,1],cb_sp,x,y,deep_face,0.,[data.vel_x,data.vel_y],data.m_p,data.h_p,L,mask,corners_face,width,r"$\log_{10}v$ (km/s)",0.,3.,cmap,sph_frame.vec2dslice)
    sph_frame.makesph_plot(fig,sp[irow,2],cb_sp,rad2d,z,deep_side,0.,[data.vel2d,data.vel_z],data.m_p,data.h_p,L,mask,corners_side,width,r"$\log_{10}v$ (km/s)",0.,3.,cmap,sph_frame.vec2dslice)

    sp[irow,2].yaxis.set_label_position("right")
    sp[irow,2].set_ylabel(labels[irow],size='x-large')


for isp in range(1,ncuts):
    sp[isp,0].remove()

cb_sp.yaxis.tick_left()
cb_sp.yaxis.set_label_position("left")

sp[ncuts-1,1].set_xlabel("pc")
sp[ncuts-1,1].set_ylabel("pc")

fig.subplots_adjust(hspace=0., wspace=0.) 
fig.tight_layout(pad=1.0,w_pad=0.0,h_pad=0.)
plt.savefig("../../figures/tcuts_montage_" + run_id + output_dir + ".png", dpi=150)
