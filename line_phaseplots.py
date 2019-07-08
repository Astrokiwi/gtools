print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os

import phaseplots

import gizmo_tools

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P

import copy

# run_name = "a2_e01"
# run_id = "3001"
# L = 300
# run_name = "newflow_vesc_thin_45"
run_name = "longrun_medflow_vesc_defaultaniso_polar"
run_id = "3032"
snap_str = "100"
L = 100

# for paper
# lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
# lineLabels = [r"CO(3-2), $866.727$ $\mu$m","CO(6-5), $433.438$ $\mu$m","HCN(4-3), $845.428$ $\mu$m","HCN(8-7), $422.796$ $\mu$m",r"H$_2$ (1-0) S(1), $2.121$ $\mu$m",r"H$_2$ (0-0) S(0), $28.18$ $\mu$m",r"H$_2$ (0-0) S(3), $9.66$ $\mu$m"]
# nx = 4
# ny = 2
# fixplort = True

#for presentation
lines = ["hcn2","co2","h2_1"]
lineLabels = [r"HCN(8-7), $422.796$ $\mu$m",r"CO(6-5), $433.438$ $\mu$m",r"H$_2$ (1-0) S(1), $2.121$ $\mu$m"]
nx = 3
ny = 1
fixplort = False

nlines = len(lines)

fig,sp = P.subplots(ny,nx,sharex=True,sharey=True,figsize=(3*nx,3*ny))
sp=sp.reshape((ny,nx))

loadVals = ["nH_p","dustTemp"]+lines

rcut = 80.
time,values = phaseplots.loadvalues(run_id,run_name,snap_str,loadVals,rcut)

for iline,line in enumerate(lines):
    ix = iline%nx
    iy = iline//nx
    this_sp = sp[iy,ix]
    mappablePlot = phaseplots.plot_phaseplot(this_sp,values,"nH_p","dustTemp",rcut,weight=line,vrange=[-20.,1.7],bins=(L,L),cmap="inferno") #,cmap='Greys'
    this_sp.set_title(lineLabels[iline])

for ix in range(nx):
    for iy in range(ny):
        iline = ix+iy*nx
        this_sp = sp[iy,ix]
        if iline>=nlines:
            this_sp.set_axis_off()
            sp[iy,ix].xaxis.set_visible(True)
        if ix!=0:
            this_sp.yaxis.set_visible(False)
        if iy!=ny-1:
            this_sp.xaxis.set_visible(False)

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
P.savefig("../figures/line_phases_EWASS2019_{}.png".format(run_name),dpi=200)
