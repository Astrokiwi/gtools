print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os

import phaseplots

import gizmo_tools

import matplotlib.pyplot as P

run_names = ["q2redo","q2edd10_aniso1redo","q2edd05redo","q2edd10redo","q2edd20redo","q2edd10_aniso3redo"]

dump_times = [100,200,500,1000]

run_id = "2014"

plotVals = ["rad_p","vel"]

rcut = 80.

if __name__ == '__main__':
    nruns = len(run_names)
    sp_lx = int(np.ceil(np.sqrt(nruns)))
    sp_ly = nruns//sp_lx

    movieDir = gizmo_tools.getMovieDir()
    for dump_time in dump_times:
        snap_str = str(dump_time)
        figtitle = "../figures/phasegrid{}{}{}{}.png".format(run_id,plotVals[0],plotVals[1],snap_str)
        
        fig,sp = P.subplots(sp_ly,sp_lx,sharex=True,sharey=True,figsize=(12,6))
        
        for run_index,run_name in enumerate(run_names):
        
            time,values = phaseplots.loadvalues(run_id,run_name,snap_str,plotVals,rcut)
        
            sp_y = run_index%sp_ly
            sp_x = run_index//sp_ly
            this_sp = sp[sp_y,sp_x]
        
            mappablePlot = phaseplots.plot_phaseplot(this_sp,values,plotVals[0],plotVals[1],rcut)
#             this_sp.set_title(run_name+"{},{}".format(sp_x,sp_y))
            if sp_x!=0:
#                 this_sp.yaxis.set_visible(False)
                this_sp.set_ylabel("")
            if sp_y!=sp_ly-1:
                this_sp.set_xlabel("")
#                 this_sp.xaxis.set_visible(False)
#             fig.colorbar(mappablePlot,label=r"$\log M$ (M$_\odot$)")


        fig.suptitle(r"$t="+("%.4f" % time)+"$ Myr", y=1.00)
        fig.tight_layout()
        fig.savefig(figtitle)
        P.close()