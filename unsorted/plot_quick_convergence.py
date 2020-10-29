import numpy as np
import gizmo_tools
import matplotlib.pyplot as plt

run_id="3032"

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
# ny = 4


titlegroups_unfilted = [  ["longrun_weakflow_settled_defaultaniso","longrun_weakflow_vesc_defaultaniso","longrun_weakflow_rapid_defaultaniso"],
                ["longrun_weakflow_settled_defaultaniso_polar","longrun_weakflow_vesc_defaultaniso_polar","longrun_weakflow_rapid_defaultaniso_polar"],
                ["longrun_medflow_settled_defaultaniso_polar","longrun_medflow_vesc_defaultaniso_polar","longrun_medflow_vesc_defaultaniso"],
                ["newflow_settled_thin_up","newflow_vesc_thin_side","newflow_vesc_thin_45","newflow_vesc_thin_up"]
                ]


# col_indices = [1,2,3,4,5,6]
col_indices = [1,4]

if __name__=='__main__':
    labels = ["coolr","coolvc","coolvr","warmr","warmvc","warmvr"]

    nx = len(col_indices)
    nruns = len(run_names)
    ny = len(titlegroups_unfilted)
    
    plot_indices = [[run_names.index(run_name) for run_name in sublist] for sublist in titlegroups_unfilted]
    
#     plot_indices = np.array_split(np.arange(nruns),ny)

#     fig,sp = plt.subplots(ny,nx,figsize=(nx*3,ny*3),sharex=True,sharey='col')
    fig,sp = plt.subplots(ny,nx,figsize=(nx*6,ny*3),sharex=True,sharey='col')

    run_data = [np.loadtxt(f"data/components_summary_{run_id}_{run_name}.dat")
                for run_name in run_names]

    for iy in range(ny):
        for ix,col_index in enumerate(col_indices):
            for plot_index in plot_indices[iy]:
#                 sp[iy,ix].plot(run_data[plot_index][:,0],run_data[plot_index][:,ix+1],label=run_names[plot_index])
                sp[iy,ix].plot(run_data[plot_index][:,0],run_data[plot_index][:,col_index],label=run_names[plot_index])
#                 sp[iy,ix].set_yscale('log')
        sp[iy,0].legend()

    for ix in range(nx):
        sp[0,ix].set_title(labels[ix])
        
    fig.savefig("../figures/convergence_check.pdf")
    

    plt.close('all')