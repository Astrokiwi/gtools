import gizmo_tools
import pynbody
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


time,snap = gizmo_tools.load_gizmo_nbody("3032","longrun_medflow_vesc_defaultaniso_polar","200")

def rot_los(snap,rot_deg):
    rot_angle = np.radians(rot_deg)
    sin_los = np.sin(rot_angle)
    cos_los = np.cos(rot_angle)

    snap["vlos"] = snap["vx"]*sin_los + snap["vy"]*cos_los

    snap["xlos"] = snap["x"]
    snap["ylos"] = snap["y"]*sin_los + snap["z"]*cos_los


snap["r2d"] = np.sqrt(snap["x"]**2+snap["y"]**2)

snap["vc"] = (snap["vx"]*snap["y"]-snap["vy"]*snap["x"])/snap["r2d"]

rot_los(snap,10.)

vert_slice = snap[np.abs(snap["xlos"])<snap["smooth"]]
horiz_slice = snap[np.abs(snap["ylos"])<snap["smooth"]]

def weighted_hist2d(snap,xval,yval,zval,**kwargs):
    MH,xedges,yedges = np.histogram2d(snap[xval],snap[yval],weights=snap["mass"],**kwargs)
    WH,xedges,yedges = np.histogram2d(snap[xval],snap[yval],weights=snap[zval],**kwargs)
    return WH/MH,xedges,yedges


fig,sp = plt.subplots(2,4,figsize=(16,8),constrained_layout=True)

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

for iay,xbound in enumerate([20.e-3,1.e-3]):
    for iax,ykey in enumerate(["vlos","vr","vc"]):
        ax = sp[iay,iax]
        hist_2d_default(ax,vert_slice["ylos"],vert_slice[ykey],'ylos',ykey,xbound,200.)
#         ax.scatter(vert_slice["ylos"],vert_slice[ykey],marker='x')
#         ax.set_xlim([-xbound,xbound])
#         ax.set_ylim([-200.,200.])
#         ax.set_xlabel('ylos')
#         ax.set_ylabel(ykey)

inner_slice = vert_slice[vert_slice["ylos"]<1.e-3]

hist_2d_default(sp[0,3],inner_slice["vc"],inner_slice["vlos"],'vc inner 1 pc','vlos inner 1 pc',200.,200.,square=True)
hist_2d_default(sp[1,3],inner_slice["vr"],inner_slice["vlos"],'vr inner 1 pc','vlos inner 1 pc',200.,200.,square=True)

fig.savefig("../figures/v_decompose.png")

plt.close('all')