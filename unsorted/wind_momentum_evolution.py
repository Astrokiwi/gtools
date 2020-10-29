import matplotlib.pyplot as plt
import numpy as np
import gizmo_tools
from scipy import stats
import sys

#     runs = [  ["longrun_weakflow_settled_defaultaniso","longrun_weakflow_vesc_defaultaniso","longrun_weakflow_rapid_defaultaniso"],
#                 ["longrun_weakflow_settled_defaultaniso_polar","longrun_weakflow_vesc_defaultaniso_polar","longrun_weakflow_rapid_defaultaniso_polar"],
#                 ["longrun_medflow_settled_defaultaniso_polar","longrun_medflow_vesc_defaultaniso","longrun_medflow_vesc_defaultaniso_polar"],
#                 ["newflow_settled_thin_up","newflow_vesc_thin_side","newflow_vesc_thin_45","newflow_vesc_thin_up"]
#                 ]

run_id="3032"
# runs=[ "longrun_weakflow_settled_defaultaniso","longrun_weakflow_vesc_defaultaniso","longrun_weakflow_rapid_defaultaniso",
#             "longrun_medflow_settled_defaultaniso_polar","longrun_medflow_vesc_defaultaniso","longrun_medflow_vesc_defaultaniso_polar"]
# run_names = runs

runs=["longrun_medflow_settled_defaultaniso_polar","longrun_medflow_vesc_defaultaniso"]
run_names = ["settled","vesc"]


# snap_strs=["000","001","005","010","050","100"]
# snap_strs=["010","100"]
snap_strs=["100"]

phi_bin_edges = np.linspace(0.,90.,25)
phi_centres = (phi_bin_edges[:-1]+phi_bin_edges[1:])/2.
phi_step = phi_bin_edges[1]-phi_bin_edges[0]

r_bin_edges = np.arange(.5,16,1.)
r_bin_centres = (r_bin_edges[:-1]+r_bin_edges[1:])/2.

# r_to_plot = [1.,5.,10.,15.]
r_to_plot = [1.,15.]

i_to_plot = np.searchsorted(r_bin_edges,r_to_plot)-1

n_rbins = len(r_to_plot)
n_tbins = len(snap_strs)

scale = 1.
fig,sp = plt.subplots(n_tbins,n_rbins,figsize=(4.*n_rbins*scale,3.*n_tbins*scale),sharex=True,sharey=True,squeeze=False)

nx = n_rbins
ny = n_tbins

times = []

for i_tbin,snap_str in enumerate(snap_strs):
    for i_run,run in enumerate(runs):
        header,data=gizmo_tools.load_gizmo_pandas(run_id,run,snap_str,["Masses","Coordinates","Velocities"])
        t = header['time']
        print("loaded time=",t)
        if i_run==0:
            times.append(t)
            print(times)
        gizmo_tools.calculate_phi(data)
        gizmo_tools.calculate_vrad(data)
        data["rad_mom"] = data["vrad"]*data["Masses"]

        binstats = stats.binned_statistic_2d(data["phi"], data["rad3d"], data["rad_mom"], statistic='sum',bins=[phi_bin_edges,r_bin_edges])
#         binstats = stats.binned_statistic_2d(data["phi"], data["rad3d"], data["Masses"], statistic='sum',bins=[phi_bin_edges,r_bin_edges])
        bin_map = binstats[0]
        bin_map/=phi_step
        for i_rbin,bin_index in enumerate(i_to_plot):
            p=sp[i_tbin,i_rbin].plot(phi_centres,bin_map[:,bin_index],label=run_names[i_run])

            outflow_slice = bin_map[:,bin_index]>0.
            outflow_slice[0] = False # don't include disk
            if np.any(outflow_slice):
                mean_phi_out = np.sum(bin_map[outflow_slice,bin_index]*phi_centres[outflow_slice])/np.sum(bin_map[outflow_slice,bin_index])
                sp[i_tbin,i_rbin].axvline(x=mean_phi_out,c=p[0].get_color()+"aa",ls='--')

sp[(len(i_to_plot)-1)//2,0].set_ylabel(r"$\mathrm{d}p/\mathrm{d}\phi$ (M$_\odot$km/s $/^\circ$)")
# sp[(len(i_to_plot)-1)//2,0].set_ylabel(r"$\mathrm{d}m/\mathrm{d}\phi$ (M$_\odot/^\circ$)")

for ix in range(nx):
    sp[ny-1,ix].set_xlabel(r"$\phi$ ($^\circ$)")
    sp[0,ix].set_title(r"${:2.1f}\mathrm{{pc}}<=r<{:2.1f}\mathrm{{pc}}$".format(r_bin_edges[i_to_plot[ix]],r_bin_edges[i_to_plot[ix]+1]))
    for iy in range(ny):
        sp[iy,ix].set_yscale('symlog', linthreshy=25.)
#         sp[iy,ix].set_yscale('log')
        sp[iy,ix].axhline(y=0.,c='k',ls='--')
        sp[iy,ix].grid(which='both',linewidth=.5,color='grey')
        sp[iy,ix].set_xlim([0.,90.])
        sp[iy,ix].set_xticks([0.,15.,30.,45.,60.,75.,90.])
sp[-1,-1].legend(loc='best',fontsize='xx-small')
        

for iy in range(ny):
    sp[iy,nx-1].yaxis.set_label_position('right')
    sp[iy,nx-1].set_ylabel(r"$t={:.4f}$ Myr".format(times[iy]))
fig.tight_layout()

# fig.savefig("../figures/wind_momentum_evolution.png",dpi=150)
fig.savefig("../figures/wind_momentum_evolution_ewass.png",dpi=150)