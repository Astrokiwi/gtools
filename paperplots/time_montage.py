print("Importing")

# Import libraries to do our magic

from ..tools import sph_frame

from sys import path
path.append("../visualisation/")
from ..tools import gizmo_tools

import matplotlib.pyplot as P

print("Running")

# output_dir = "q2redo"
# runs = ["q2edd05redo","q2edd10_aniso1redo","q2edd10_aniso3redo","q2edd10redo","q2edd20redo","q2redo"]
# runs = ["q2redo"]
# run_id = "2014"
# 
# snapx = [0,100,200,500,1000]
# snapx = [200,500]


# runs = ["run_a2_e01"]
# run_id = "2022"

runs = ["a2_e01"]
run_id = "3001"
snapx = [0,10,20,50,100]



nsnaps = len(snapx)

gizmoDir = gizmo_tools.getGizmoDir(run_id)
movieDir = gizmo_tools.getMovieDir()

flattenedPlot = True
rot = [0.,0.]
plot_thing = ['temp','nH']
L=256
width = .8
corners_face = [-width/2.,-width/2.]
corners_side = [0.,-width/2.]

temp_cmap = "plasma"
nH_cmap = "viridis"

for output_dir in runs:
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    infiles = [fullDir+"/snapshot_"+("%03d" % x)+".hdf5" for x in snapx]
    
    infile = infiles[3]
    ntimes = len(infiles)
    
    plotLabel, plotRanges, plotSliceTypes, plotCustomMass, plotData, logSliceTypes, extraBarTypes, plusMinusTypes, divergingTypes, customCmaps, customCmaps2 = sph_frame.pack_dicts()


    fig,sp = P.subplots(ntimes,6,figsize=(12.,12.), gridspec_kw = {'width_ratios':[1,16,16,16,16,1],'height_ratios':[16]*ntimes})
    cb_sp = sp[0,0]
    cb_sp_temp = sp[0,5]
#     fig,sp = P.subplots(ntimes+1,2,figsize=(6.,3.*ntimes), gridspec_kw = {'height_ratios':[16]*ntimes+[1]})
#     cb_sp = sp[ntimes,0]

    for irow in range(ntimes-1):
        for icol in [1,2,3,4]:
            sp[irow,icol].set_xticklabels([])
            sp[irow,icol].set_xticklabels([])
            sp[irow,icol].get_shared_x_axes().join(sp[irow,icol], sp[irow+1,icol])

        for icol in [1,2]:
            sp[irow,icol].get_shared_y_axes().join(sp[irow,icol], sp[irow,icol+2])    

    for irow in range(ntimes):
        for icol in [2,3,4]:
            sp[irow,icol].set_yticklabels([])
        for icol in [1,2,3,4]:
            sp[irow,icol].set(adjustable='box-forced', aspect='equal')
   
    for irow,infile in enumerate(infiles):
        time,data,x,y,z,rad2d,deep_face,deep_side,mask,n,n_ones = sph_frame.load_process_gadget_data(infile,rot,plot_thing,plotData,ringPlot=flattenedPlot,flatPlot=flattenedPlot)
        sph_frame.makesph_plot(fig,sp[irow,1],cb_sp,x,y,deep_face,0.,data.nH_p,data.m_p,data.h_p,L,mask,corners_face,width,r"$\log_{10} n_{H}$ (cm$^{-3}$)",0.,8.,nH_cmap,sph_frame.weightslice,contour=False)
        sph_frame.makesph_plot(fig,sp[irow,2],cb_sp,rad2d,z,deep_side,0.,data.nH_p,data.m_p,data.h_p,L,mask,corners_side,width,r"$\log_{10} n_{H}$ (cm$^{-3}$)",0.,8.,nH_cmap,sph_frame.weightslice,contour=False)
        sph_frame.makesph_plot(fig,sp[irow,3],cb_sp_temp,x,y,deep_face,0.,data.TK_p,data.m_p,data.h_p,L,mask,corners_face,width,r"$\log_{10} T_g$ (K)",0.,5.,temp_cmap,sph_frame.weightslice,contour=False)
        sph_frame.makesph_plot(fig,sp[irow,4],cb_sp_temp,rad2d,z,deep_side,0.,data.TK_p,data.m_p,data.h_p,L,mask,corners_side,width,r"$\log_{10} T_g$ (K)",0.,5.,temp_cmap,sph_frame.weightslice,contour=False)
        
        sp[irow,4].yaxis.set_label_position("right")
        sp[irow,4].set_ylabel("t={0:.2f} kyr".format(time*1.e3),size='x-large')

    for isp in range(1,ntimes):
        sp[isp,0].remove()
        sp[isp,5].remove()
#     sp[ntimes,1].remove()    
#     P.colorbar(sp[0,0].collections[0],ax=cb_sp,orientation='horizontal')
    cb_sp.yaxis.tick_left()
    cb_sp.yaxis.set_label_position("left")
    
    sp[ntimes-1,1].set_xlabel("pc")
    sp[ntimes-1,1].set_ylabel("pc")
    
    fig.subplots_adjust(hspace=0., wspace=0.) 
    fig.tight_layout(pad=0.0,w_pad=0.0,h_pad=0.)
    P.savefig("../../figures/time_montage_"+run_id+output_dir+".png",dpi=150)
    #     outfile = "../../figures/nHevolve"+run_id+output_dir+".png"
    #     sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['nH']*nsnaps,L=400,scale=4.)
    # 
    #     outfile = "../../figures/Tevolve"+run_id+output_dir+".png"
    #     sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['temp']*nsnaps,L=400,scale=4.)
