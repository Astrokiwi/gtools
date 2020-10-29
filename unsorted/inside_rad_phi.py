import numpy as np
import matplotlib.pyplot as plt

from unsorted import parameter_mangler
import gizmo_tools
    
if __name__=='__main__':
    run_id = "2030"

#     snap_str = "010"
#     for snap_str in "001","004","010"
#     snap_str = "001"
        
    manglia = parameter_mangler.load_manglia("prams/pram_options_testflows.dat")
    key_combinations = parameter_mangler.mangle_combinations(manglia)
#     print(key_combinations)
    run_dirs = [parameter_mangler.filebase_from_combination("testflows", x) for x in key_combinations]
#     print(run_dirs)
    key_bases = [list(x.mangle_texts.keys()) for x in manglia]
    mangle_n = [len(x) for x in key_bases]
    n_mangles = len(mangle_n)
    
    iaxis = 1
    nx = mangle_n[iaxis]
    jaxis = 0
    ny = mangle_n[jaxis]
    color_axis = 2
    colors = ['r','g','b','k']
#     linestyles = ['dashed','dashdot','solid']

    for snap_str in "001","004","010":
        fig,sp = plt.subplots(ny,nx,figsize=(4*nx,4*ny),sharex=True,sharey=True,squeeze=False)
        phi_bin_edges = np.linspace(0.,90.,19)
        phi_centres = (phi_bin_edges[:-1]+phi_bin_edges[1:])/2.
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
                t,data=gizmo_tools.load_gizmo_pandas(run_id,run_dir,snap_str,["Masses","Coordinates"])
                print(run_dir,"found")
                gizmo_tools.calculate_phi(data)
                phi_indices = np.digitize(data["phi"],phi_bin_edges)-1
                phi_slices = [(phi_indices==i) for i in range(len(phi_bin_edges)-1)]
                nbin = np.array([np.sum(slice) for slice in phi_slices])
                data["rad3d"] = np.sqrt(data["Coordinates_x"]**2+data["Coordinates_y"]**2+data["Coordinates_z"]**2)
                n_inside = np.array([np.sum(data["rad3d"][slice]<1.) for slice in phi_slices])
                inside_frac = n_inside/nbin
                print(ix,iy,sp.shape)
                sp[ix,iy].plot(phi_centres,1.-inside_frac,label=color_label,color=colors[ic])
            except OSError:
                print(run_dir," not found")
                continue
        x = np.linspace(0.,90.,100)
        y = np.zeros(x.shape)
        for ix in range(nx):
            for iy in range(ny):
                sp[ix,iy].set_title(key_bases[jaxis][iy]+"_"+key_bases[iaxis][ix])
        sp[0,0].legend(loc='best',fontsize='xx-small')
        for sp_row in sp:
            for p in sp_row:
                p.legend(loc='best',fontsize='xx-small')
#                 p.plot(x,y,linestyle=':')
                p.set_xticks(phi_bin_edges)
                p.set_ylim(-.1,1.1)
#                 p.set_yscale("symlog",linthreshy=100.)
#                 p.set_yticks(ytick_values)
#                 p.set_yticklabels(ytick_values)
        plt.savefig("../figures/inside_rad_phi_testflows{}.png".format(snap_str),dpi=150)
        plt.close()
        