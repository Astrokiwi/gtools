import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

import parameter_mangler
import gizmo_tools
from scipy import stats
    
if __name__=='__main__':
    run_id = "2030"

    snap_str = "010"
        
    manglia = parameter_mangler.load_manglia("prams/pram_options_testflows.dat")
    key_combinations = parameter_mangler.mangle_combinations(manglia)
#     print(key_combinations)
    run_dirs = [parameter_mangler.filebase_from_combination("testflows",x) for x in key_combinations]
#     print(run_dirs)
    key_bases = [list(x.mangle_texts.keys()) for x in manglia]
    mangle_n = [len(x) for x in key_bases]
    n_mangles = len(mangle_n)
    

    phi_bin_edges = np.linspace(0.,90.,19)
    phi_centres = (phi_bin_edges[:-1]+phi_bin_edges[1:])/2.
    
    r_bin_edges = np.linspace(.5,7.5,8)
    r_bin_centres = (r_bin_edges[:-1]+r_bin_edges[1:])/2.
    n_r_bins = r_bin_centres.size

    iaxis = 1
    nx = mangle_n[iaxis]
    jaxis = 0
    ny = mangle_n[jaxis]
    color_axis = 2
    colors = ['r','g','b','k']
#     linestyles = ['dashed','dashdot','solid']
    
    fig = np.empty((n_r_bins),dtype=object)
    sp = np.empty((n_r_bins,nx,ny),dtype=object)
    
    for ifig in range(n_r_bins):
         fig[ifig],sp[ifig,:,:] = plt.subplots(ny,nx,figsize=(4*nx,4*ny),sharex=True,sharey=True,squeeze=False)
#     for irun,run_dir in enumerate(run_dirs):
    for irun,run_dir in enumerate(run_dirs[:34]):
        try:
            coords = np.zeros((n_mangles),dtype=int)
            for ii in range(n_mangles):
                if ii==n_mangles-1:
                    coords[ii] = irun%mangle_n[ii]
                coords[ii] = (irun//np.prod(mangle_n[ii+1:]))%mangle_n[ii]
            ix = coords[iaxis]
            iy = coords[jaxis]
            ic = coords[color_axis]
            color_label = key_bases[color_axis][ic]
            t,data=gizmo_tools.load_gizmo_pandas(run_id,run_dir,snap_str,["Masses","Coordinates","Velocities"])
            print(run_dir,"found")
            gizmo_tools.calculate_phi(data)
            gizmo_tools.calculate_vrad(data)
            data["rad3d"] = np.sqrt(data["Coordinates_x"]**2+data["Coordinates_y"]**2+data["Coordinates_z"]**2)
#             binstats = stats.binned_statistic_2d(data["phi"], data["rad3d"], data["vrad"], statistic='mean',bins=[phi_bin_edges,r_bin_edges])
            binstats = stats.binned_statistic_2d(data["phi"], data["rad3d"], data["vrad"]*data["Masses"], statistic='sum',bins=[phi_bin_edges,r_bin_edges])
            bin_map = binstats[0]
            print(bin_map.shape)
            for ifig in range(n_r_bins):
                print(ifig,ix,iy,sp.shape)
                sp[ifig,ix,iy].plot(phi_centres,bin_map[:,ifig],label=color_label,color=colors[ic])
        except OSError as e:
            print(run_dir," not found")
            continue
    x = np.linspace(0.,90.,100)
    y = np.zeros(x.shape)
#     ytick_values = [-50,-20,0,20,50,70,100,200,500,700]
    for ifig in range(n_r_bins):
        for ix in range(nx):
            for iy in range(ny):
                sp[ifig,ix,iy].set_title(key_bases[jaxis][iy]+"_"+key_bases[iaxis][ix])
        sp[ifig,0,0].legend(loc='best',fontsize='xx-small')
        for sp_row in sp[ifig,:,:]:
            for p in sp_row:
#                 p.legend(loc='best',fontsize='xx-small')
                p.plot(x,y,linestyle=':')
                p.set_xticks(phi_bin_edges[::2])
#                 p.set_ylim(ytick_values[0],ytick_values[-1])
                p.set_ylim(-1.e5,1.e5)
                p.set_yscale("symlog",linthreshy=100.)
#                 p.set_yscale("symlog",linthreshy=10.)
#                 p.set_yscale("log")
#                 p.set_yscale("symlog",linthreshy=50.)
#                 p.set_yticks(ytick_values)
#                 p.set_yticklabels(ytick_values)
        fig[ifig].suptitle(r"$%g\mathrm{pc}<=r<%g\mathrm{pc}$" % (r_bin_edges[ifig],r_bin_edges[ifig+1]) )
        fig[ifig].text(0.5, 0.04, r"$\phi$ ($^\circ$)", ha='center')
#         fig[ifig].text(0.04, 0.5, r"$v_r$ (km s$^{-1}$)", va='center', rotation='vertical')
        fig[ifig].text(0.04, 0.5, r"$\Sigma m v_r$ (M$_\odot$ km s$^{-1}$)", va='center', rotation='vertical')
#         fig[ifig].savefig("../figures/vrad_rad_phi_testflows_r{}.png".format(ifig),dpi=150)
        fig[ifig].savefig("../figures/prad_rad_phi_testflows_r{}.png".format(ifig),dpi=150)
    plt.close('all')
        