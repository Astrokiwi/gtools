print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os

import phaseplots

import gizmo_tools

import matplotlib.pyplot as P


# run_ids = ["2014"]
# run_names = ["q2redo"]
# snap_strs = ["1000"]

# run_ids = ["2014","2014","2014","1055","1055","2018"]
# run_names = ["q2redo","q2redo","q2edd10redo","q2_SN","q2_SN_slow","treetest"]
# # snap_strs = ["511","1000","511","511","511"]
# snap_strs = ["132","1000","132","132","132","132"]
# titles = ["Run A","Run A","Run C2",r"Run A$_{**}$",r"Run A$_{*}$","treefix"]

# doruns = [True,False,True,True,False]
# doruns = [True,False,True,True,True,True]

# run_names = ["run_a0_e1","run_a1_e1","run_a2_e01","run_a2_e01_T1000","run_a2_e01_T300","run_a2_e01_T30","run_a2_e05","run_a2_e1","run_a2_e2","run_a3_e1"]
# run_names = ["run_a2_e01"]
run_names = ["a2_e01"]
# run_ids = ["2022"]*len(run_names)
run_ids = ["3001"]*len(run_names)
snap_strs = ["020"]*len(run_names)
doruns = [True]*len(run_names)
# titles = [name[4:] for name in run_names]
titles = run_names

# run_names = ["newtable00001","newtable00005","newtable","newtable0005","newtable001","newtable005","newtable01"]
# run_ids = ["2022"]*len(run_names)
# snap_strs = ["039"]*len(run_names)
# doruns = [True]*len(run_names)
# titles = run_names

for run_id,run_name,snap_str,title,dorun in zip(run_ids,run_names,snap_strs,titles,doruns):

    if not dorun:
        continue

    loadVals = ["dustTemp","vradpoly","TK_p"]
#     loadVals = ["TK_p","vradpoly","TK_p"]
    coords = [ [0,0], [1,0], [1,1] ]
    plotPairs = [ [0,1], [0,2], [1,2] ]


    rcut = 80.

    fig,sp = P.subplots(2,2,sharex='col',sharey='row',figsize=(6,7))

    time,values = phaseplots.loadvalues(run_id,run_name,snap_str,loadVals,rcut)

    print("time=",time)

    for iplot,plotVals in enumerate(plotPairs):
        this_sp=sp[coords[iplot][0],coords[iplot][1]]
        mappablePlot = phaseplots.plot_phaseplot(this_sp,values,loadVals[plotVals[0]],loadVals[plotVals[1]],rcut,cmap='Greys')

    sp[0,1].set_axis_off()
    # sp[0,2].set_axis_off()

    sp[0,0].axvline(x=2.,c='k',ls='--')
    sp[1,0].axvline(x=2.,c='k',ls='--')
    sp[1,0].axhline(y=3.,c='k',ls='--')
    sp[1,1].axhline(y=3.,c='k',ls='--')
    
    

    sp[0,0].set_xlabel("")
    
    P.suptitle(title+" t=%3.1f kyr"%(time*1.e3))
    sp[1,1].set_ylabel("")

    # fig.subplots_adjust(right=0.80)
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    fig.subplots_adjust(bottom=0.20,left=.2,right=.95,top=.9)
    cbar_ax = fig.add_axes([0.15, 0.07, 0.7, 0.03])

    cbar = P.colorbar(mappablePlot,label=r"$\log M$ (M$_\odot$)",cax=cbar_ax,orientation='horizontal')

    # fig.tight_layout()
    P.savefig("../figures/phaseplot_vrad_temp{}{}{}.png".format(run_id,run_name,snap_str),dpi=200)
