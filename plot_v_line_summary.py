print("setting up")

import numpy as np
import matplotlib.pyplot as plt
import gizmo_tools
from sys import exit


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

# run_names = ["longrun_medflow_vesc_defaultaniso_polar"]

nruns = len(run_names)

lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3","dust"]
line_order = [0,1,2,3,4,6,5,7]
# nlines = len(line_order)

brightness_keys = ["lum_"+line for line in lines]

line_labels = gizmo_tools.lineNamesFull.copy()
line_labels['dust'] = "$F_{IR}$"

print("loading")
v_data = dict()
for run_name in run_names:
    v_data[run_name] = np.loadtxt(f"data/v_summary_{run_name}.dat")

print("Setting up parameters")
run_parameters = gizmo_tools.load_run_parameters("3032")

run_parameters = {x:y for x,y in run_parameters.items() if x in run_names}

ordered_keys = gizmo_tools.run_parameters_names(run_parameters)
# ordered_keys = [x for x in ordered_keys_in if x in run_names]

if False:
    print("plotting lines for each run")
    #first pass - plot all lines for each run
    nx,ny = gizmo_tools.gridsize_from_n(nruns)
    fig,sp = plt.subplots(ny,nx,figsize=(4.*nx,4.*ny),sharex=True,sharey=True,squeeze=False,constrained_layout=True)
    for irun,run_name in enumerate(ordered_keys):
        ix = irun%nx
        iy = irun//nx
        ax = sp[iy,ix]
    
        ax.set_title(run_parameters[run_name]['name'])
        ax.set_yscale('symlog',linthreshy=1.e-2)

        for iline,line_name in enumerate(brightness_keys):
            v_circ = v_data[run_name][iline,0]
            v_rad = v_data[run_name][iline,1]
            ax.scatter(v_circ,v_rad,label=line_name,edgecolor=gizmo_tools.colors[irun%gizmo_tools.ncolors],marker=gizmo_tools.markers[iline//gizmo_tools.ncolors])

    for ix in range(nx):
        ax = sp[ny-1,ix]
        ax.set_xlabel(r"$v_c$ (km/s)")
    for iy in range(ny):
        ax = sp[iy,0]
        ax.set_ylabel(r"$v_r$ (km/s)")
    #     ax.set_ylabel(r"$\log_{10}r_{90}$ (pc)")
    sp[0,0].legend()
    for iplot in range(len(run_names),nx*ny):
        ix = iplot%nx
        iy = iplot//nx
        ax = sp[iy,ix]
        ax.set_axis_off()
    fig.savefig("../figures/line_v_lines.pdf")

    plt.close('all')







print("plotting runs for each line")
#second pass - plot all runs for each line

print("vc vs vr")
nlums = len(brightness_keys)
nx,ny = gizmo_tools.gridsize_from_n(nlums,aspect=2)
line_rows = ny
ny*=3
scale = .7
fig,sp = plt.subplots(ny,nx,figsize=(3.5*nx*scale,3.*ny*scale),sharex=True,sharey=True,squeeze=False,constrained_layout=True)
for ext_offset,ext_suffix in enumerate(["","ext_face_","ext_edge_"]):
#     fig,sp = plt.subplots(ny,nx,figsize=(3.*nx,3.*ny),sharex=True,sharey=True,squeeze=False,constrained_layout=True)
    for iline in range(nlums):
#     for iline,line_name in enumerate(lines):
        ix = iline%nx
        iy = iline//nx
        ax = sp[iy+line_rows*ext_offset,ix]
    
        datarow = line_order[iline] + ext_offset*nlums
        line_name = lines[line_order[iline]]

    #     label_key = line_name[4:]
        ax.set_title(line_labels[line_name])
        ax.set_yscale('symlog',linthreshy=1.e0)
        ax.set_xscale('symlog',linthreshy=1.e0)

        x = []
        y = []

#         full_legend_plot = ix==0 and iy==0 and ext_offset==0
        full_legend_plot = ix==3 and iy==1 and ext_offset==2

        for irun,run_name in enumerate(ordered_keys):
            v_circ = v_data[run_name][datarow,0]
            v_rad = v_data[run_name][datarow,1]
            label = run_parameters[run_name]['name'] if full_legend_plot else None
            ax.scatter(v_circ,v_rad,label=label,color=run_parameters[run_name]['color'],marker=run_parameters[run_name]['marker'],s=run_parameters[run_name]['size'],linewidth=run_parameters[run_name]['thickness'])
            x.append(v_circ)
            y.append(v_rad)
        fit = np.polyfit(x,y,1)
        x = np.linspace(0.,np.max(x),100)
        y = fit[0]*x+fit[1]
        ax.plot(x,y,ls='--',label=r"$v_r={:.2f}{:+.2f}v_c$".format(fit[1],fit[0]))
    #     print(fit)
#         ax.legend(loc='best',fontsize=4 if full_legend_plot else 'xx-small',ncol=2 if full_legend_plot else 1)
        ax.legend(loc='best',fontsize='xx-small' if full_legend_plot else 'small',ncol=2 if full_legend_plot else 1)

for ix in range(nx):
    ax = sp[ny-1,ix]
    ax.set_xlabel(r"$v_c$ (km/s)")
for iy in range(ny):
    ax = sp[iy,0]
    ax.set_ylabel(r"$v_r$ (km/s)")
    # for iplot in range(len(brightness_keys)):
    #     ix = iplot%nx
    #     iy = iplot//nx
    #     ax = sp[iy,ix]
    #     x = np.linspace(0.,50.,100)
    #     y = 500.-x*10.
    #     ax.plot(x,y,ls='--',label=r"$v_r=B-Av_x$ fit")

# for iplot in range(len(brightness_keys),nx*ny):
#     ix = iplot%nx
#     iy = iplot//nx
#     ax = sp[iy,ix]
#     ax.set_axis_off()
fig.savefig(f"../figures/line_vr_vc_runs.pdf")
plt.close('all')






if False:
    print("edge-on vdisp vs vcirc")
    nlums = len(brightness_keys)
    nx,ny = gizmo_tools.gridsize_from_n(nlums,aspect=2)
    line_rows = ny
    ny*=2
    fig,sp = plt.subplots(ny,nx,figsize=(3.*nx,3.*ny),sharex=True,sharey=True,squeeze=False,constrained_layout=True)
    for ext_offset,ext_suffix in enumerate(["","ext_edge_"]):
        for iline,line_name in enumerate(lines):
            ix = iline%nx
            iy = iline//nx
            ax = sp[iy+line_rows*ext_offset,ix]
    
            datarow = iline + ext_offset*2*nlums

        #     label_key = line_name[4:]
            ax.set_title(line_labels[line_name])
            ax.set_yscale('symlog',linthreshy=1.e0)
            ax.set_xscale('symlog',linthreshy=1.e0)

            x = []
            y = []

            for irun,run_name in enumerate(ordered_keys):
                v_circ = v_data[run_name][datarow,0]
                disp_y = v_data[run_name][datarow,2]
                label = run_parameters[run_name]['name'] if ix==0 and iy==0 else None
                ax.scatter(v_circ,disp_y,label=label,color=run_parameters[run_name]['color'],marker=run_parameters[run_name]['marker'],s=run_parameters[run_name]['size'],linewidth=run_parameters[run_name]['thickness'])
                x.append(v_circ)
                y.append(disp_y)
            fit = np.polyfit(x,y,1)
            x = np.linspace(0.,np.max(x),100)
            y = fit[0]*x+fit[1]
            ax.plot(x,y,ls='--',label=r"$\sigma_y={:.2f}{:+.2f}v_c$".format(fit[1],fit[0]))
        #     print(fit)
            ax.legend(loc='best',fontsize='xx-small')

    for ix in range(nx):
        ax = sp[ny-1,ix]
        ax.set_xlabel(r"$v_c$ (km/s)")
    for iy in range(ny):
        ax = sp[iy,0]
        ax.set_ylabel(r"$\sigma_y$ (km/s)")
        # for iplot in range(len(brightness_keys)):
        #     ix = iplot%nx
        #     iy = iplot//nx
        #     ax = sp[iy,ix]
        #     x = np.linspace(0.,50.,100)
        #     y = 500.-x*10.
        #     ax.plot(x,y,ls='--',label=r"$v_r=B-Av_x$ fit")

    # for iplot in range(len(brightness_keys),nx*ny):
    #     ix = iplot%nx
    #     iy = iplot//nx
    #     ax = sp[iy,ix]
    #     ax.set_axis_off()
    fig.savefig(f"../figures/line_vdisp_vc_runs.pdf")
    plt.close('all')
