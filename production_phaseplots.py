print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os

import phaseplots

import gizmo_tools

import matplotlib.pyplot as P


run_id = "2014"

run_name = "q2redo"

snap_str = "1000"

loadVals = ["dustTemp","vrad","TK_p"]
coords = [ [0,0], [1,0], [1,1] ]
plotPairs = [ [0,1], [0,2], [1,2] ]


rcut = 80.

fig,sp = P.subplots(2,2,sharex='col',sharey='row',figsize=(8,7))

time,values = phaseplots.loadvalues(run_id,run_name,snap_str,loadVals,rcut)

for iplot,plotVals in enumerate(plotPairs):
    this_sp=sp[coords[iplot][0],coords[iplot][1]]
    mappablePlot = phaseplots.plot_phaseplot(this_sp,values,loadVals[plotVals[0]],loadVals[plotVals[1]],rcut,cmap='Greys')

sp[0,1].set_axis_off()
# sp[0,2].set_axis_off()

sp[0,0].axvline(x=30.,c='k',ls='--')
sp[1,0].axvline(x=30.,c='k',ls='--')
sp[1,0].axhline(y=10**2.5,c='k',ls='--')
sp[1,1].axhline(y=10**2.5,c='k',ls='--')

sp[0,0].set_xlabel("")
sp[1,1].set_ylabel("")

fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])

cbar = P.colorbar(mappablePlot,label=r"$\log M$ (M$_\odot$)",cax=cbar_ax)

# fig.tight_layout()
P.savefig("../figures/phaseplot_vrad_temp.png",dpi=200)
