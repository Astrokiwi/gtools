print("Importing")

# Import libraries to do our magic
import numpy as np
import sys

from tools import phaseplots

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P

# run_name = "a2_e01"
# run_id = "3001"
# L = 300
# run_name = "newflow_vesc_thin_45"
run_name = "longrun_medflow_vesc_defaultaniso_polar"
run_id = "3032"
snap_str = "200"
# L = 100
L = 64

# for paper
# lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
# lineLabels = [r"CO(3-2), $866.727$ $\mu$m","CO(6-5), $433.438$ $\mu$m","HCN(4-3), $845.428$ $\mu$m","HCN(8-7), $422.796$ $\mu$m",r"H$_2$ (1-0) S(1), $2.121$ $\mu$m",r"H$_2$ (0-0) S(0), $28.18$ $\mu$m",r"H$_2$ (0-0) S(3), $9.66$ $\mu$m"]
# nx = 4
# ny = 2
# fixplort = True

#for presentation
# lines = ["hcn2","co2","h2_1"]
# lineLabels = [r"HCN(8-7), $422.796$ $\mu$m",r"CO(6-5), $433.438$ $\mu$m",r"H$_2$ (1-0) S(1), $2.121$ $\mu$m"]
# nx = 3
# ny = 1
# fixplort = False

# for big paper, for information
lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
# lines = ["co1","h2_1"]
lineLabels = [r"CO(3-2), $866.727$ $\mu$m","CO(6-5), $433.438$ $\mu$m","HCN(4-3), $845.428$ $\mu$m","HCN(8-7), $422.796$ $\mu$m",r"H$_2$ (1-0) S(1), $2.121$ $\mu$m",r"H$_2$ (0-0) S(0), $28.18$ $\mu$m",r"H$_2$ (0-0) S(3), $9.66$ $\mu$m"]
nx = 4
ny = 2
fixplort = False


nlines = len(lines)

loadVals = ["nH_p","dustTemp","TK_p"]+lines

# rcut = 80.
rcut = 800.
time,values = phaseplots.loadvalues(run_id, run_name, snap_str, loadVals, rcut)


print("one axis charts")

n_dependence = True

if n_dependence:
#     ncuts = 3
    ncuts = 3
else:
    ncuts = 1

nH_cut = 1.e3
tdust_cut=15.
tdust_cut_big=100.


# line plots
fig,sp = P.subplots(nlines,3*ncuts,sharey=True,sharex='col',figsize=(3*3*ncuts,3*nlines))
for icut in range(ncuts):
    for iline,line in enumerate(lines):
        for ival,val in enumerate(["nH_p","dustTemp","TK_p"]):
            this_sp = sp[iline,ival+icut*3]
            if n_dependence:
                if icut==0:
                    val_slice = values["dustTemp"]<tdust_cut
                    title=r"$T_d<15$ K"
                elif icut==1:
                    val_slice = (values["dustTemp"]>=tdust_cut) & (values["dustTemp"]<tdust_cut_big)
                    title=r"$15\leq T_d<100$ K"
                else:
                    val_slice = (values["dustTemp"]>=tdust_cut_big)
                    title=r"$T_d\geq100$ K"

# 
#                 if ncuts==2:
#                     if icut==0:
#                         val_slice = values["dustTemp"]<tdust_cut
#                         title=r"$T_d<15$ K"
#                     else:
#                         val_slice = values["dustTemp"]>=tdust_cut
#                         title=r"$T_d\geq15$ K"
#                 elif ncuts==3:
#                     if icut==0:
#                         val_slice = values["nH_p"]<nH_cut
#                         title=r"$nH_p<10^2$ cm$^{-2}$"
#                     elif icut==1:
#                         val_slice = (values["nH_p"]>=nH_cut) & (values["dustTemp"]<tdust_cut)
#                         title=r"$nH_p\geq10^2$ cm$^{-2}$, $T_d<10^2$ K"
#                     else:
#                         val_slice = (values["nH_p"]>=nH_cut) & (values["dustTemp"]>=tdust_cut)
#                         title=r"$nH_p\geq10^2$ cm$^{-2}$, $T_d\geq10^2$ K"
                if iline==0:
                    this_sp.set_title(title)
            else:
                val_slice = np.ones(len(values["nH_p"]),dtype=np.bool)
            mappablePlot,x,y = phaseplots.plot_phaseplot(this_sp, values, val, line, rcut, bins=(L, L), cmap="inferno", slice=val_slice, drawmed=True)
            outfile_name = f"data/median_emmissivity_{run_id}_{run_name}_{snap_str}_{val}_{line}_{icut}.dat"
            print(outfile_name)
            np.savetxt(outfile_name,np.array([x,y]).T)

        
#         arg_order = np.argsort(values[val])
#         sp[ival].loglog(values[val][arg_order],values[line][arg_order],label=lineLabels[iline])
#     sp[1].loglog(values["dustTemp"],values[line])
#     sp[2].loglog(values["TK_p"],values[line])
# sp[0].legend()
P.savefig("../figures/line_phases_chart_{}.png".format(run_name))
P.close('all')

sys.exit()


print("phase plots")
fig,sp = P.subplots(ny*2,nx,sharex=True,sharey='row',figsize=(3*nx,3*ny*2),squeeze=False)
# fig,sp = P.subplots(ny,nx,sharex=True,sharey=True,figsize=(3*nx,3*ny))
# sp=sp.reshape((ny,nx))
for iline,line in enumerate(lines):
    ix = iline%nx
    iy = iline//nx
    this_sp = sp[iy,ix]
    mappablePlot = phaseplots.plot_phaseplot(this_sp,values,"nH_p","dustTemp",rcut,weight=line,vrange=[-20.,1.7],bins=(L,L),cmap="inferno") #,cmap='Greys'
    this_sp.set_title(lineLabels[iline])

    iy+=ny
    this_sp = sp[iy,ix]
    mappablePlot = phaseplots.plot_phaseplot(this_sp,values,"nH_p","TK_p",rcut,weight=line,vrange=[-20.,1.7],bins=(L,L),cmap="inferno") #,cmap='Greys'
    this_sp.set_title(lineLabels[iline])


# for ix in range(nx):
#     for iy in range(ny):
#         iline = ix+iy*nx
#         this_sp = sp[iy,ix]
#         if iline>=nlines:
#             this_sp.set_axis_off()
#             sp[iy,ix].xaxis.set_visible(True)
#             sp[iy+ny,ix].xaxis.set_visible(True)
#         if ix!=0:
#             this_sp.yaxis.set_visible(False)
#             sp[iy+ny,ix].yaxis.set_visible(False)
#         if iy!=ny-1:
#             this_sp.xaxis.set_visible(False)
#             sp[iy+ny,ix].xaxis.set_visible(False)

P.colorbar(mappablePlot,label="Line Emissivity (erg/s/g)")

if fixplort:
    cbax = fig.axes[-1]
    cbax_bbox = cbax.get_position()
    cbar_width = cbax_bbox.x1-cbax_bbox.x0
    cbar_height = cbax_bbox.y1-cbax_bbox.y0
    old_bbox = fig.axes[-2].get_position() # get colorbar location
    cbax.set_position([old_bbox.x0,old_bbox.y0,cbar_width,cbar_height])
else:
    fig.tight_layout()
# P.savefig("../figures/line_phases.png",dpi=200)
# P.savefig("../figures/line_phases_TORUS2018.png",dpi=200)
# P.savefig("../figures/line_phases_EWASS2019_{}.png".format(run_name),dpi=200)
P.savefig("../figures/line_phases_big_{}.png".format(run_name),dpi=200)
P.close('all')



