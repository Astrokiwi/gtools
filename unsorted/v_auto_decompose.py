import gizmo_tools
import pynbody
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches


zoom_boxes = [
                 [[-1.5,1.5],[50.,150.]]
                ,[[-1.5,1.5],[-150.,-50.]]
                
            ]


time,snap = gizmo_tools.load_gizmo_nbody("3032","longrun_medflow_vesc_defaultaniso_polar","200")

def rot_los(snap,rot_deg):
    rot_angle = np.radians(rot_deg)
    sin_los = np.sin(rot_angle)
    cos_los = np.cos(rot_angle)

    snap["vlos"] = snap["vx"]*sin_los + snap["vy"]*cos_los

    snap["xlos"] = snap["x"]
    snap["ylos"] = snap["y"]*sin_los + snap["z"]*cos_los

    snap["zlos"] = snap["x"]*sin_los + snap["y"]*cos_los


snap["r2d"] = np.sqrt(snap["x"]**2+snap["y"]**2)

snap["vc"] = (snap["vx"]*snap["y"]-snap["vy"]*snap["x"])/snap["r2d"]

rot_los(snap,10.)

vert_slice = snap[np.abs(snap["xlos"])<snap["smooth"]]
horiz_slice = snap[np.abs(snap["ylos"])<snap["smooth"]]

def weighted_hist2d(snap,xval,yval,zval,**kwargs):
    MH,xedges,yedges = np.histogram2d(snap[xval],snap[yval],weights=snap["mass"],**kwargs)
    WH,xedges,yedges = np.histogram2d(snap[xval],snap[yval],weights=snap[zval],**kwargs)
    return WH/MH,xedges,yedges

ny = len(zoom_boxes)

for slice,slice_axis in [(vert_slice,"ylos"),(horiz_slice,"xlos")]:

    fig,sp = plt.subplots(ny,4,figsize=(16,ny*4),constrained_layout=True,squeeze=False)

    def hist_2d_default(ax,x,y,xlabel,ylabel,xbound,ybound,square=False):
        H,xedges,yedges = np.histogram2d(x,y,range=[[-xbound,xbound],[-ybound,ybound]],bins=(100,100))
        ax.pcolormesh(xedges,yedges,H.T,norm=colors.LogNorm())
    #     ax.scatter(x,y,marker='x')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    #     if square:
    #         ax.set_xlim(-xbound,xbound)
    #         ax.set_ylim(-ybound,ybound)
    #         ax.axis('equal')

    print(slice_axis,len(slice))

    # slice=vert_slice
    # slice_axis = "ylos"

    for iy,zoombox in enumerate(zoom_boxes):

        zoom_box_kpc = np.array(zoombox[0])*1.e-3
        zoom_box_vlos = np.array(zoombox[1])

        ax=sp[iy,0]
        hist_2d_default(ax,slice[slice_axis],slice["vlos"],slice_axis,"vlos",20.e-3,200.)
        rect = patches.Rectangle((zoom_box_kpc[0],zoom_box_vlos[0]),zoom_box_kpc[1]-zoom_box_kpc[0],zoom_box_vlos[1]-zoom_box_vlos[0],linewidth=1,edgecolor='r',facecolor='none')
        ax.add_patch(rect)
    
        box_snap = slice[(slice[slice_axis]>=zoom_box_kpc[0]) & (slice[slice_axis]<=zoom_box_kpc[1]) & (slice["vlos"]>=zoom_box_vlos[0]) & (slice["vlos"]<zoom_box_vlos[1])]
    
        print(len(box_snap))
        
        for iplot,(yval,ybound) in enumerate([
                             ["vlos",200.]
                            ,["vr",200.]
                            ,["vc",200.]
                            ]):
            ax = sp[iy,iplot+1]
            hist_2d_default(ax,box_snap["zlos"],box_snap[yval],r"$z_{los}$ (kpc)",yval,5.e-3,ybound)
    #         hist_2d_default(ax,box_snap["y"],box_snap[yval],r"$z_{derot}$ (kpc)",yval,5.e-3,ybound)




    fig.savefig("../figures/v_auto_decompose_{}.png".format(slice_axis))

    plt.close('all')