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

#https://gist.github.com/thriveth/8560036 - colourblind colours from here
colors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']
ncolors=len(colors)
cmap = mpl.cm.get_cmap('coolwarm')

markers = ["x","+",'|']

rad_data = dict()
for run_name in run_names:
    rad_data[run_name] = np.loadtxt("data/summary_lumrads_{}.dat".format(run_name))

run_parameters = gizmo_tools.load_run_parameters("3032")

run_parameters = {x:y for x,y in run_parameters.items() if x in run_names}

ordered_keys = gizmo_tools.run_parameters_names(run_parameters)


#first pass - plot all lines for each run
nx,ny = gizmo_tools.gridsize_from_n(nruns)
fig,sp = plt.subplots(ny,nx,figsize=(4.*nx,4.*ny),sharex=True,sharey=True,squeeze=False,constrained_layout=True)
for irun,run_name in enumerate(ordered_keys):
    ix = irun%nx
    iy = irun//nx
    ax = sp[iy,ix]
    
    ax.set_title(run_parameters[run_name]['name'])

    for iline,line_name in enumerate(brightness_keys):
        r_50 = rad_data[run_name][iline,0]
        r_90 = rad_data[run_name][iline,1]
        log_r_50 = np.log10(r_50)
        log_r_90_50 = np.log10(r_90/r_50)
        log_r_90 = np.log10(r_90)
        ax.scatter(log_r_50,log_r_90_50,label=line_name,edgecolor=gizmo_tools.colors[irun%ncolors],marker=gizmo_tools.markers[iline//ncolors])
#         ax.annotate(line_name,(log_r_90_50,log_r_50))
#         ax.scatter(log_r_90_50,log_r_90,label=line_name)
#         ax.annotate(line_name,(log_r_90_50,log_r_90))
# for ix in range(nx):
#     for iy in range(ny):
#         ax = sp[iy,ix]
#         ax.set_xscale('log')
#         ax.set_yscale('log')
for ix in range(nx):
    ax = sp[ny-1,ix]
    ax.set_xlabel(r"$\log_{10}r_{50}$ (pc)")
for iy in range(ny):
    ax = sp[iy,0]
    ax.set_ylabel(r"$\log_{10}r_{90}/r_{50}$")
#     ax.set_ylabel(r"$\log_{10}r_{90}$ (pc)")
sp[0,0].legend()
for iplot in range(len(run_names),nx*ny):
    ix = iplot%nx
    iy = iplot//nx
    ax = sp[iy,ix]
    ax.set_axis_off()

fig.savefig("../figures/line_r_50_90_lines.pdf")

plt.close('all')

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

    for irun,run_name in enumerate(ordered_keys):
        r_50 = rad_data[run_name][iline,0]
        r_90 = rad_data[run_name][iline,1]
        log_r_50 = np.log10(r_50)
        log_r_90 = np.log10(r_90)
        log_r_90_50 = np.log10(r_90/r_50)
#         ax.scatter(log_r_50,log_r_90_50,label=run_name,edgecolor=colors[irun%ncolors],marker=markers[irun//ncolors])
#         ax.scatter(log_r_50,log_r_90_50,label=run_name,color=run_parameters[run_name]['color'],marker=run_parameters[run_name]['marker'],s=run_parameters[run_name]['size'],linewidth=run_parameters[run_name]['thickness'])
        ax.scatter(log_r_50,log_r_90,label=run_parameters[run_name]['name'],color=run_parameters[run_name]['color'],marker=run_parameters[run_name]['marker'],s=run_parameters[run_name]['size'],linewidth=run_parameters[run_name]['thickness'])
#         ax.scatter(log_r_90_50,log_r_90,label=run_name)
#         ax.annotate(run_name,(log_r_50,log_r_90_50),rotation=270.)
for ix in range(nx):
    ax = sp[ny-1,ix]
    ax.set_xlabel(r"$\log_{10}r_{50}$ (pc)")
for iy in range(ny):
    ax = sp[iy,0]
#     ax.set_ylabel(r"$\log_{10}r_{90}/r_{50}$")
    ax.set_ylabel(r"$\log_{10}r_{90}$ (pc)")
for iplot in range(len(brightness_keys)):
    ix = iplot%nx
    iy = iplot//nx
    ax = sp[iy,ix]
    x = np.linspace(0.,2.5,100)
    #         y = 1.5-x
    #         ax.plot(x,y,ls='--',label=r"Constant $r_{90}$")

    y = x
    ax.plot(x,y,ls='--',label=r"$r_{90}\propto r_{50}$")

for iplot in range(len(brightness_keys),nx*ny):
    ix = iplot%nx
    iy = iplot//nx
    ax = sp[iy,ix]
    ax.set_axis_off()


sp[0,0].legend(loc='best',fontsize='xx-small')
# plt.legend()
#         ax.set_ylim([-0.1,None])

#     ax.set_ylabel(r"$\log_{10}r_{90}$ (pc)")
# for ix in range(nx):
#     for iy in range(ny):
#         ax = sp[iy,ix]
#         ax.set_xscale('log')
#         ax.set_yscale('log')
fig.savefig("../figures/line_r_50_90_runs.pdf")

plt.close('all')
