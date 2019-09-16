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

nruns = len(run_names)

lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]

brightness_keys = ["luminosity"]+["lum_"+line for line in lines]

line_labels = gizmo_tools.lineNamesFull.copy()
line_labels['nosity'] = "$F_{IR}$"

print("loading")
v_data = dict()
for run_name in run_names:
    v_data[run_name] = np.loadtxt(f"data/v_summary_{run_name}.dat")

print("Setting up parameters")
run_parameters = gizmo_tools.load_run_parameters("3032")

run_parameters = {x:y for x,y in run_parameters.items() if x in run_names}

ordered_keys = gizmo_tools.run_parameters_names(run_parameters)
# ordered_keys = [x for x in ordered_keys_in if x in run_names]

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
nlums = len(brightness_keys)
nx,ny = gizmo_tools.gridsize_from_n(nlums,aspect=2)
fig,sp = plt.subplots(ny,nx,figsize=(3.*nx,3.*ny),sharex=True,sharey=True,squeeze=False,constrained_layout=True)
for iline,line_name in enumerate(brightness_keys):
    ix = iline%nx
    iy = iline//nx
    ax = sp[iy,ix]
    
    label_key = line_name[4:]
    ax.set_title(line_labels[label_key])
    ax.set_yscale('symlog',linthreshy=1.e0)

    x = []
    y = []

    for irun,run_name in enumerate(ordered_keys):
        v_circ = v_data[run_name][iline,0]
        v_rad = v_data[run_name][iline,1]
        label = run_parameters[run_name]['name'] if ix==0 and iy==0 else None
        ax.scatter(v_circ,v_rad,label=label,color=run_parameters[run_name]['color'],marker=run_parameters[run_name]['marker'],s=run_parameters[run_name]['size'],linewidth=run_parameters[run_name]['thickness'])
        x.append(v_circ)
        y.append(v_rad)
    fit = np.polyfit(x,y,1)
    x = np.linspace(0.,np.max(x),100)
    y = fit[0]*x+fit[1]
    ax.plot(x,y,ls='--',label=r"$v_r={:.2f}{:+.2f}v_c$".format(fit[1],fit[0]))
#     print(fit)
    ax.legend(loc='best',fontsize='xx-small')
    
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

for iplot in range(len(brightness_keys),nx*ny):
    ix = iplot%nx
    iy = iplot//nx
    ax = sp[iy,ix]
    ax.set_axis_off()
fig.savefig("../figures/line_v_runs.pdf")

plt.close('all')
