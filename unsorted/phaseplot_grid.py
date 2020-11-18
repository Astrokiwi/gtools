print("Importing")

# Import libraries to do our magic
import numpy as np

from tools import phaseplots

import gizmo_tools

import matplotlib.pyplot as P

# run_names = ["q2redo","q2edd10_aniso1redo","q2edd05redo","q2edd10redo","q2edd20redo","q2edd10_aniso3redo"]
# dump_times = [100,200,500,1000]
# run_id = "2014"

# run_names = ["q2redo","q2edd10_aniso1redo","q2edd05redo","q2edd10redo","q2edd20redo","q2edd10_aniso3redo"]
# dump_times = [10,20,50,100]
# run_id = "2022"

# run_names = ["SN_test_m1","SN_test_m2","SN_test_m3","SN_test_m3_thinner","SN_test_m3_thin","SN_test_T5_half","SN_test_T5","SN_test_T6_half","SN_test_T6","SN_test_T7_half","SN_test_T7"]
# dump_times = [10,20,30]
# run_id = "2030"


# plot_title = "SF_comparison"
# run_names_2030 = ["SN_test_m2_slowlong","SN_test_m3_slowlong","SN_test_m2_slowerlong","SN_test_m3_slowerlong"]
# run_names_3030 = ["no_SNe_long"]
# run_names = run_names_2030 + run_names_3030
# dump_times = [10,20,30]
# run_ids = ["2030"]*len(run_names_2030)+["3030"]*len(run_names_3030)
# 
# 
# # plotVals = ["rad_p","vel"]
# plotVals = ["nH_p","TK_p"]

plot_title = "testflows"
run_names = ["testflows_medflow_settled_defaultaniso","testflows_medflow_settled_doubleiso","testflows_medflow_settled_iso",
             "testflows_medflow_settled_thiniso","testflows_medflow_vesc_defaultaniso","testflows_medflow_vesc_doubleiso",
             "testflows_medflow_vesc_thiniso","testflows_weakflow_rapid_defaultaniso","testflows_weakflow_rapid_doubleiso",
             "testflows_weakflow_rapid_iso","testflows_weakflow_rapid_thiniso","testflows_weakflow_settled_defaultaniso",
             "testflows_weakflow_settled_doubleiso","testflows_weakflow_settled_iso","testflows_weakflow_settled_thiniso",
             "testflows_weakflow_vesc_defaultaniso","testflows_weakflow_vesc_doubleiso","testflows_weakflow_vesc_iso",
             "testflows_weakflow_vesc_thiniso"]

dump_times = [10]
run_ids = ["2030"]*len(run_names)


# plotVals = ["rad_p","vel"]
plotVals = ["nH_p","TK_p","tau"]

rcut = 20.

if __name__ == '__main__':
    nruns = len(run_names)
    sp_lx = int(np.ceil(np.sqrt(nruns)))
    sp_ly = int(np.ceil(nruns/sp_lx))
    
    xPlotVal = plotVals[0]

    movieDir = gizmo_tools.getMovieDir()
    for dump_time in dump_times:
        for yPlotVal in plotVals[1:]:
            snap_str = "%03d"%dump_time
            figtitle = "../figures/phasegrid{}{}{}{}.png".format(plot_title,xPlotVal,yPlotVal,snap_str)
        
            fig,sp = P.subplots(sp_ly,sp_lx,sharex=True,sharey=True,figsize=(sp_lx*4.,sp_ly*4.))
            title_time = -1
        
            for run_index,(run_id,run_name) in enumerate(zip(run_ids,run_names)):
                try:
                    time,values = phaseplots.loadvalues(run_id, run_name, snap_str, [xPlotVal, yPlotVal], rcut)
                except OSError:
                    print("directory"+run_id+" simulation "+run_name+" dump "+snap_str+" not found, skipping")
                    continue # skip this file
        
                title_time = time
            
                sp_y = run_index%sp_ly
                sp_x = run_index//sp_ly
                this_sp = sp[sp_y,sp_x]
        
                mappablePlot = phaseplots.plot_phaseplot(this_sp, values, xPlotVal, yPlotVal, rcut)
    #             this_sp.set_title(run_name+"{},{}".format(sp_x,sp_y))
                this_sp.set_title(run_name)
                if sp_x!=0:
    #                 this_sp.yaxis.set_visible(False)
                    this_sp.set_ylabel("")
                if sp_y!=sp_ly-1:
                    this_sp.set_xlabel("")
    #                 this_sp.xaxis.set_visible(False)
    #             fig.colorbar(mappablePlot,label=r"$\log M$ (M$_\odot$)")


            fig.suptitle(r"$t="+("%.4f" % title_time)+"$ Myr", y=1.00)
            fig.tight_layout()
            fig.savefig(figtitle)
            P.close()