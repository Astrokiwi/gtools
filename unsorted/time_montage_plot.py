print("Importing")

# Import libraries to do our magic
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pickle
import gizmo_tools

print("Running")

test_mode = False

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

run_parameters = gizmo_tools.load_run_parameters("3032")
# run_parameters = {x:y for x,y in run_parameters.items() if x in runs}
gizmo_tools.run_parameters_names(run_parameters)

run_id = "3032"
snapx = [10,20,50,100,150,200,250,300]

chosen_snaps = [2,3,4,5,6,7]

nruns = len(runs)
ntimes = len(chosen_snaps)

flattenedPlot = True
rot = [0.,0.]
plot_thing = ['temp','nH','vels']
L=256
# L=32
width = 20.
corners_face = [-width/2.,-width/2.]
corners_side = [0.,-width/2.]

xedges = np.linspace(0.,width,L+1)
yedges = np.linspace(-width/2.,width/2.,L+1)

xmids = (xedges[:-1]+xedges[1:])/2.
ymids = (yedges[:-1]+yedges[1:])/2.

temp_cmap = "plasma"
nH_cmap = "viridis"
vel_cmap = "plasma"

# scale = 1.
scale = 0.8
fig_titles = ["nH_side","Tg_side","v_side"]
nfig = len(fig_titles)
# figures = [plt.subplots(ntimes,nruns+1,figsize=((3.*nruns+1.)*scale,(3.*ntimes)*scale), gridspec_kw = {'width_ratios':[16]*nruns+[1],'height_ratios':[16]*ntimes}, sharex='col',sharey='col',constrained_layout=True) for x in range(nfig)]
# fig_titles = ["nH_face","nH_side","Tg_face","Tg_side"]

# figs = [plt.figure(figsize=((3.*nruns+1.)*scale,(3.*ntimes)*scale),constrained_layout=True) for x in range(nfig)]
figs = [plt.figure(figsize=((3.*nruns+1.)*scale,(3.*ntimes)*scale)) for x in range(nfig)]
gridspecs = [gridspec.GridSpec(ntimes,nruns*16+2,figure=fig) for fig in figs]
subplots = []
for fig,gs in zip(figs,gridspecs):
    sp = np.empty((ntimes,nruns+1),dtype=np.object)
    for itime in range(ntimes):
        for irun in range(nruns):
#             share_ax_x = None if itime==0 else sp[0,irun]
#             share_ax_y = None if irun==0 else sp[itime,0]
            share_ax_x = None #sp[0,0]
            share_ax_y = None #sp[0,0]
            sp[itime,irun] = fig.add_subplot(gs[itime:(itime+1),irun*16:(irun+1)*16],sharex=share_ax_x,sharey=share_ax_y)
        sp[itime,-1]=fig.add_subplot(gs[itime:(itime+1),nruns*16+1:nruns*16+2])
    subplots.append(sp)
figures = [[fig,sp] for fig,sp in zip(figs,subplots)]

cb_ix = nruns

# set up plots
for fig,sp in figures:
#     for irow in range(ntimes-1):
#         for icol in range(nruns):
#             sp[irow,icol].set_xticklabels([])
#             sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow+1,icol])
#             sp[irow,icol].get_shared_y_axes().join(sp[irow,icol], sp[irow+1,icol])
# 
#         for icol in range(nruns-1):
#             sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow,icol+1])
#             sp[irow,icol].get_shared_y_axes().join(sp[irow,icol], sp[irow,icol+1])
            
    for irow in range(ntimes):
        for icol in range(1,nruns):
            sp[irow,icol].set_yticklabels([])
    for irow in range(ntimes):
        for icol in range(nruns):
            sp[irow,icol].set_aspect(aspect='equal')
#             sp[irow,icol].set_aspect(aspect='equal',adjustable='box-forced')
    sp[ntimes-1,0].set_xlabel("pc")
    sp[ntimes-1,0].set_ylabel("pc")

def load_maybe(filename):
    if test_mode:
        time = irow*1.e5
        if "_v_" in filename:
            sph_map = [np.random.random((L,L))]*3
        else:
            sph_map = np.random.random((L,L))
    else:
        map_data = pickle.load(open(filename,"rb"))
        time = map_data[0]
        sph_map = map_data[1]
    return time,sph_map

def plot_run(irun,output_dir):
    for fig,sp in figures:
        sp[0,irun].set_title(run_parameters[output_dir]["name"])#,fontsize=6)

    for irow,snap_index in enumerate(chosen_snaps):
        try:
            fig,sp = figures[0]
            time,sph_map = load_maybe(f"data/time_montage_frame_n_{run_id}_{output_dir}_{snap_index}.dat")
            cmappable = sp[irow,irun].pcolormesh(xedges,yedges,sph_map,cmap=nH_cmap,vmin=-2.,vmax=8.)
            fig.colorbar(cmappable,label=r"$\log_{10} n_{H}$ (cm$^{-3}$)",cax=sp[irow,cb_ix])

            fig,sp = figures[1]
            time,sph_map = load_maybe(f"data/time_montage_frame_T_{run_id}_{output_dir}_{snap_index}.dat")
            cmappable = sp[irow,irun].pcolormesh(xedges,yedges,sph_map,cmap=temp_cmap,vmin=0.,vmax=5.)
            fig.colorbar(cmappable,label=r"$\log_{10} T_g$ (K)",cax=sp[irow,cb_ix])

            fig,sp = figures[2]
            time,sph_maps = load_maybe(f"data/time_montage_frame_v_{run_id}_{output_dir}_{snap_index}.dat")
            sph_map,map1,map2 = sph_maps
            cmappable = sp[irow,irun].pcolormesh(xedges,yedges,sph_map,cmap=vel_cmap,vmin=0.,vmax=3.)
            sp[irow,irun].streamplot(xmids,ymids,map1,map2)
            fig.colorbar(cmappable,label=r"$\log_{10} v$ (km/s)",cax=sp[irow,cb_ix])

    #         fig,sp = figures[2]
    #         sph_map = np.loadtxt(f"data/time_montate_frame_V_{run_id}_{output_dir}_{irow}.dat")

            for fig,sp in figures:
                sp[irow,nruns-1].yaxis.set_label_position("right")
                sp[irow,nruns-1].set_ylabel("t={:.0f} kyr".format(time*1.e3))#,size='x-large')

            if irow<ntimes-1:
                for fig,sp in figures:
                    sp[irow,irun].set_xticklabels([])
            if irun>0:
                for fig,sp in figures:
                    sp[irow,irun].set_yticklabels([])

        except OSError as e:
            print(output_dir,irow," can't be opened, leaving panel blank")
            for fig,sp in figures:
                sp[irow,irun].set_visible(False)
#             return
    
def add_zoom():
    fig,sp = figures[0]
#     time,sph_map = pickle.load(open(f"data/time_montage_zoom.dat","rb"))
#     pos1 = sp[-1,1].get_position()
#     pos2 = sp[-1,2].get_position()
#     x0 = (pos1.x0+pos2.x0)/2
#     x1 = (pos1.x1+pos2.x1)/2
#     y0 = pos1.y0
#     y1 = pos1.y1
#     ax = fig.add_axes([x0,y0,x1-x0,y1-y0])
    ax = fig.add_subplot(gs[ntimes-1:ntimes,20:36])
    cax = fig.add_subplot(gs[ntimes-1:ntimes,37:38])
    time,sph_map = load_maybe(f"data/time_montage_zoom.dat")
#     cmappable = ax.pcolormesh(np.linspace(0.,10.,256+1),np.linspace(0.,10.,256+1),sph_map,cmap=nH_cmap,vmin=0.,vmax=5.)
    cmappable = ax.pcolormesh(np.linspace(0.,10.,512+1),np.linspace(0.,10.,512+1),sph_map,cmap=nH_cmap,vmin=0.,vmax=5.)
    ax.set_aspect(aspect='equal')
    fig.colorbar(cmappable,label=r"$\log_{10} n_{H}$ (cm$^{-3}$)",cax=cax)


    gizmo_tools.box_connected_two_axes(sp[3,1]
                        ,ax
                        ,[.1,.1]
                        ,[9.8,9.8]
                        ,color='red')
#                         ,color='#10D7AE')


if __name__ == '__main__':
    # do plots
    for irun,output_dir in enumerate(runs):
        plot_run(irun,output_dir)
    add_zoom()

    

    #dirty hack to add zoom
#     fig,sp = figures[0]
#     map_data = pickle.load(open(f"data/time_montage_zoom.dat","rb"))
#     sph_map = map_data[1]
#     cmappable = sp[irow,irun].pcolormesh(np.linspace(0.,10.,512+1),np.linspace(0.,10.,512+1),sph_map,cmap=nH_cmap,vmin=0.,vmax=5.)
#     fig.colorbar(cmappable,label=r"$\log_{10} n_{H}$ (cm$^{-3}$)",cax=sp[irow,cb_ix])

    for fig_title,(fig,sp) in zip(fig_titles,figures):
        fig.subplots_adjust(hspace=0.1, wspace=0.0,left=0.05,right=0.95,top=0.95,bottom=0.05)
#         fig.subplots_adjust(hspace=0., wspace=0.)
#         fig.tight_layout(pad=0.0,w_pad=0.0,h_pad=0.)
#         fig.tight_layout()
        fig.savefig("../figures/time_montage_together_simple_"+fig_title+"_"+run_id+".png",dpi=150)
