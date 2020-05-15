print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
from multiprocessing import Pool

#from sys import path
#path.append("../src/")
import sph_frame

from sys import path
path.append("../")
import gizmo_tools

import matplotlib.pyplot as plt

from functools import partial

import pickle

print("Running")

run_id = "3032"
# runs = [
# "longrun_medflow_settled_defaultaniso",
# "longrun_medflow_settled_defaultaniso_polar",
# "longrun_medflow_vesc_defaultaniso",
# "longrun_medflow_vesc_defaultaniso_polar",
# "longrun_weakflow_rapid_defaultaniso",
# "longrun_weakflow_rapid_defaultaniso_polar",
# "longrun_weakflow_settled_defaultaniso",
# "longrun_weakflow_settled_defaultaniso_polar",
# "longrun_weakflow_vesc_defaultaniso",
# "longrun_weakflow_vesc_defaultaniso_polar",
# "newflow_settled_thin_up",
# "newflow_vesc_thin_45",
# "newflow_vesc_thin_side",
# "newflow_vesc_thin_up"]
runs = [
"longrun_weakflow_vesc_defaultaniso_polar",
"longrun_medflow_vesc_defaultaniso_polar",
"newflow_vesc_thin_45",
"newflow_vesc_thin_up"
]
snapx = [10,20,50,100,150,200,250,300]
width = 20.


# snapx = [10,20,50,100,150,200]


run_id = "3032"
# snapx = [0,10,20,50,100,150,200,250,300]
snapx = [100]
# run_id = "3032"
# runs = [
# "longrun_medflow_vesc_defaultaniso_polar",
# "longrun_medflow_vesc_defaultaniso_polar_highedd"]
# #"longrun_medflow_vesc_defaultaniso_polar_highedd_rsub"]
# snapx = [200]
# width = 40.



nruns = len(runs)
ntimes = len(snapx)
gizmoDir = gizmo_tools.getGizmoDir(run_id)
movieDir = gizmo_tools.getMovieDir()

flattenedPlot = True
rot = [0.,0.]
plot_thing = ['temp','nH','vels']
L=256
# L=32
# width = 200.
# corners_face = [-width/2.,-width/2.]

ringPlot=True
corners_side = [0.,-width/2.]

# ringPlot=False
# corners_side = [-width/2.,-width/2.]

temp_cmap = "plasma"
nH_cmap = "viridis"
vel_cmap = "plasma"

scale = 1.
fig_titles = ["nH_side","Tg_side","v_side"]
nfig = len(fig_titles)
figures = [plt.subplots(ntimes,nruns+1,figsize=((3.*nruns+1.)*scale,(3.*ntimes)*scale), gridspec_kw = {'width_ratios':[16]*nruns+[1],'height_ratios':[16]*ntimes},squeeze=False) for x in range(nfig)]
# fig_titles = ["nH_face","nH_side","Tg_face","Tg_side"]

cb_ix = nruns

# set up plots
for fig,sp in figures:
    for irow in range(ntimes-1):
        for icol in range(nruns):
            sp[irow,icol].set_xticklabels([])
            sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow+1,icol])
            sp[irow,icol].get_shared_y_axes().join(sp[irow,icol], sp[irow+1,icol])

        for icol in range(nruns-1):
            sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow,icol+1])
            sp[irow,icol].get_shared_y_axes().join(sp[irow,icol], sp[irow,icol+1])
            
    for irow in range(ntimes):
        for icol in range(1,nruns):
            sp[irow,icol].set_yticklabels([])
        for icol in range(nruns):
            sp[irow,icol].set_aspect(aspect='equal')#,adjustable='box-forced')
    sp[ntimes-1,0].set_xlabel("pc")
    sp[ntimes-1,0].set_ylabel("pc")

def plot_run(irun,output_dir):
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    infiles = [fullDir+"/snapshot_"+("%03d" % x)+".hdf5" for x in snapx]

    plotLabel, plotRanges, plotSliceTypes, plotCustomMass, plotData, logSliceTypes, extraBarTypes, plusMinusTypes, divergingTypes, symLogSliceTypes, customCmaps, customCmaps2 = sph_frame.pack_dicts()
       
    if False:
        for fig,sp in figures:
            sp[0,irun].set_title(output_dir,fontsize=6)

    for irow,infile in enumerate(infiles):
        try:
            time,data,x,y,z,rad2d,deep_face,deep_side,mask,n,n_ones = sph_frame.load_process_gadget_data(infile,rot,plot_thing,plotData,ringPlot=ringPlot,flatPlot=flattenedPlot)
        except OSError as e:
            print(infile," can't be opened, leaving panel blank")
            return
    
    #         fig,sp = figures[0]
    #         sph_frame.makesph_plot(fig,sp[irow,irun],sp[irow,cb_ix],x,y,deep_face,0.,data.nH_p,data.m_p,data.h_p,L,mask,corners_face,width,r"$\log_{10} n_{H}$ (cm$^{-3}$)",-2.,8.,nH_cmap,sph_frame.weightslice,contour=False)
    #         fig,sp = figures[2]
    #         sph_frame.makesph_plot(fig,sp[irow,irun],sp[irow,cb_ix],x,y,deep_face,0.,data.TK_p,data.m_p,data.h_p,L,mask,corners_face,width,r"$\log_{10} T_g$ (K)",0.,5.,temp_cmap,sph_frame.weightslice,contour=False)

        outmap = []

        fig,sp = figures[0]
#         sph_map = sph_frame.makesph_plot(fig,sp[irow,irun],sp[irow,cb_ix],rad2d,z,deep_side,0.,data.nH_p,data.m_p,data.h_p,L,mask,corners_side,width,r"$\log_{10} n_{H}$ (cm$^{-3}$)",-2.,8.,nH_cmap,sph_frame.weightslice,contour=False)
        sph_map = sph_frame.makesph_plot(fig,sp[irow,irun],sp[irow,cb_ix],rad2d,z,deep_side,0.,data.nH_p,data.m_p,data.h_p,L,mask,corners_side,width,r"$\log_{10} n_{H}$ (cm$^{-3}$)",-2.,6.,nH_cmap,sph_frame.weightslice,contour=False)
        pickle.dump([time,sph_map]
                    ,open(f"../data/time_montage_frame_n_{run_id}_{output_dir}_{irow}.dat","wb"))

#         np.savetxt(f"../data/time_montage_frame_n_{run_id}_{output_dir}_{irow}.dat",sph_map)

        fig,sp = figures[1]
        sph_map = sph_frame.makesph_plot(fig,sp[irow,irun],sp[irow,cb_ix],rad2d,z,deep_side,0.,data.TK_p,data.m_p,data.h_p,L,mask,corners_side,width,r"$\log_{10} T_g$ (K)",0.,5.,temp_cmap,sph_frame.weightslice,contour=False)
        pickle.dump([time,sph_map]
                    ,open(f"../data/time_montage_frame_T_{run_id}_{output_dir}_{irow}.dat","wb"))

#         np.savetxt(f"../data/time_montage_frame_T_{run_id}_{output_dir}_{irow}.dat",sph_map)

        fig,sp = figures[2]
        sph_maps = sph_frame.makesph_plot(fig,sp[irow,irun],sp[irow,cb_ix],rad2d,z,deep_side,0.,[data.vel2d,data.vel_z],data.m_p,data.h_p,L,mask,corners_side,width,r"$\log_{10} v$ (km/s)",0.,3.,vel_cmap,sph_frame.vec2dslice_flat,contour=False)
        pickle.dump([time,sph_maps]
                    ,open(f"../data/time_montage_frame_v_{run_id}_{output_dir}_{irow}.dat","wb"))

#         np.savetxt(f"../data/time_montage_frame_v_{run_id}_{output_dir}_{irow}.dat",sph_map)

        for fig,sp in figures:
            sp[irow,nruns-1].yaxis.set_label_position("right")
            sp[irow,nruns-1].set_ylabel("t={0:.2f} kyr".format(time*1.e3),size='x-large')

# with Pool(processes=176) as pool:
#     pool.starmap(plot_run
#                 ,zip(range(len(runs)),runs)
#                     )


# do plots
for irun,output_dir in enumerate(runs):
    plot_run(irun,output_dir)

for fig_title,(fig,sp) in zip(fig_titles,figures):
    fig.subplots_adjust(hspace=0., wspace=0.) 
    fig.tight_layout(pad=0.0,w_pad=0.0,h_pad=0.)
    fig.savefig("../../figures/time_montage_together_"+fig_title+"_"+run_id+".png",dpi=150)
#     fig.savefig("../../figures/time_montage_Ledd_"+fig_title+"_"+run_id+".png",dpi=300)
