import matplotlib.pyplot as plt
import numpy as np
import gizmo_tools
from scipy import stats
import sys

run_id="3032"
# run_dir="longrun_weakflow_vesc_defaultaniso_polar"
runs = [
"longrun_medflow_settled_defaultaniso",
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

# runs = runs[0:2] # slice for speed

# snaps = range(25,325,25)
snaps = range(0,325,50)
snap_strs=["{:03d}".format(x) for x in snaps]

# r_bin_edges = np.linspace(.7,1.5,5)
# r_bin_edges = np.arange(0.5,100.5,1.)                                                                                                                                                      
r_bin_edges = np.arange(0.5,50.5,1.)                                                                                                                                                      
r_bin_centres = (r_bin_edges[:-1]+r_bin_edges[1:])/2.
n_r_bins = r_bin_centres.size
dr = r_bin_edges[1]-r_bin_edges[0]

# plt.figure()
# plt.yscale('symlog',linthreshy=1.e-1)

nruns = len(runs)
ny = int(np.floor(np.sqrt(nruns)))
nx = nruns//ny
while nx*ny<nruns:
    nx+=1
scale = 1.5
fig,sp = plt.subplots(ny,nx,figsize=(4.*nx*scale,4.*ny*scale),sharex=True,sharey=True,squeeze=False)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']                                                                                                     

for irun,run in enumerate(runs):
    ix = irun%nx
    iy = irun//nx
    cur_sp = sp[iy,ix]
    
    cur_sp.set_title(run,fontsize=6.)
    cur_sp.set_yscale('symlog',linthreshy=1.e-2)
    for itime,snap_str in enumerate(snap_strs):
        try:
            header,data=gizmo_tools.load_gizmo_pandas(run_id,run,snap_str,["Masses","Coordinates","Velocities"])
        except OSError as e:
            print(snap_str," not found, skipping")
            continue
        print("{:32s}, loaded time={:4.2f} Myr".format(run,header["time"]))
        gizmo_tools.calculate_vrad(data)

        km_s_to_pc_yr = 1.02269032e-6
        outslice = (data["vrad"]>0.)
        label_string = r"$t={:4.2f}$ Myr".format(header["time"])
        if np.sum(outslice)>0:
            mdot_out_binned,bin_edges = np.histogram(data["rad3d"][outslice],weights=data["vrad"][outslice]*data["Masses"][outslice]*km_s_to_pc_yr/dr,bins=r_bin_edges)
            cur_sp.plot(r_bin_centres,mdot_out_binned,label=label_string,color=colors[itime])
            label_string = None
        if np.sum(~outslice)>0:
            mdot_in_binned,bin_edges = np.histogram(data["rad3d"][~outslice],weights=data["vrad"][~outslice]*data["Masses"][~outslice]*km_s_to_pc_yr/dr,bins=r_bin_edges)
            cur_sp.plot(r_bin_centres,mdot_in_binned,label=label_string,color=colors[itime])
        # plt.hist(bin_edges[:-1],bin_edges,weights=mdot_binned)

    cur_sp.legend(loc='best',fontsize='xx-small',ncol=3)
plt.savefig("../figures/mdot_hist_quick.png")
plt.close()