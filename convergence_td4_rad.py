import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle
from sys import exit
import gizmo_tools

run_names = [
"longrun_medflow_settled_defaultaniso_polar",
"longrun_medflow_vesc_defaultaniso",
"longrun_medflow_vesc_defaultaniso_polar",
"longrun_weakflow_rapid_defaultaniso",
"longrun_weakflow_rapid_defaultaniso_polar",
"longrun_weakflow_settled_defaultaniso",
"longrun_weakflow_settled_defaultaniso_polar",
"longrun_weakflow_vesc_defaultaniso",
"longrun_weakflow_vesc_defaultaniso_polar",
"newflow_settled_thin_up",
"newflow_vesc_thin_45",
"newflow_vesc_thin_side",
"newflow_vesc_thin_up"]

nruns = len(run_names)

lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]

brightness_keys = ["luminosity"]+["lum_"+line for line in lines]

line_labels = gizmo_tools.lineNamesFull.copy()
line_labels['nosity'] = "$F_{IR}$"

run_parameters = gizmo_tools.load_run_parameters("3032")

run_parameters = {x:y for x,y in run_parameters.items() if x in run_names}

ordered_keys = gizmo_tools.run_parameters_names(run_parameters)

scale = 2.

def line_convergences(suffix="extinguished",snap_strs=["050","075","100"]):
    rad_data = dict()
    for run_name in run_names:
        run_data = []
        for snap_str in snap_strs:
            try:
                run_data.append(np.loadtxt("data/summary_lumrads_{}_{}_{}.dat".format(run_name,suffix,snap_str)))
            except OSError:
                break # i.e. run until file not found
        rad_data[run_name] = np.array(run_data)
        print(np.array(rad_data[run_name]).shape)
    
    nx,ny = gizmo_tools.gridsize_from_n(nruns,aspect=2)
    fig,sp = plt.subplots(ny,nx,figsize=(3.*nx*scale,3.*ny*scale),sharex=True,sharey=True,squeeze=False,constrained_layout=True)

    for irun,run_name in enumerate(ordered_keys):
        ax = sp[irun//nx,irun%nx]
        x_pts = range(rad_data[run_name].shape[0])
        
        for iline,line in enumerate(brightness_keys):
            if iline==0: continue
            ax.plot(x_pts,rad_data[run_name][:,iline,0]/rad_data[run_name][8,iline,0],color=gizmo_tools.colors[iline],label=line)
            ax.plot(x_pts,rad_data[run_name][:,iline,1]/rad_data[run_name][8,iline,1],color=gizmo_tools.colors[iline],ls='--')
        
        ax.set_yscale('log')
        ax.set_title(run_name+"\n"+run_parameters[run_name]['name'])
        
    sp[0,0].legend(loc='best',fontsize='xx-small')
    
    fig.savefig("../figures/lum_rad_evolution.pdf")

if __name__=='__main__':
    line_convergences("unextinguished"
                        ,["{:03d}".format(x) for x in range(10,350,10)])
