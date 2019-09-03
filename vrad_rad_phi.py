import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

import parameter_mangler
import gizmo_tools
from scipy import stats

import h5py
    

run_id = "2030"
snap_str = "010"

mode = "v"

if __name__=='__main__':
    
    if mode!="v" and mode!="p":
        raise Exception("mode must be v (velocity) or p (momentum)")
    
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
    
    df = pd.DataFrame()
    
    df["phi"] = phi_centres
    
    valid_rundirs = []
    
    
    for irun,run_dir in enumerate(run_dirs):
#     for irun,run_dir in enumerate(run_dirs[:34]):
        try:
            t,data=gizmo_tools.load_gizmo_pandas(run_id,run_dir,snap_str,["Masses","Coordinates","Velocities"])
            print(run_dir,"found")
            gizmo_tools.calculate_phi(data)
            gizmo_tools.calculate_vrad(data)
            data["rad3d"] = np.sqrt(data["Coordinates_x"]**2+data["Coordinates_y"]**2+data["Coordinates_z"]**2)
            if mode=="v":
                binstats = stats.binned_statistic_2d(data["phi"], data["rad3d"], data["vrad"], statistic='mean',bins=[phi_bin_edges,r_bin_edges])
            elif mode=="p":
                binstats = stats.binned_statistic_2d(data["phi"], data["rad3d"], data["vrad"]*data["Masses"], statistic='sum',bins=[phi_bin_edges,r_bin_edges])
            bin_map = binstats[0]
            print(bin_map.shape)

            for ir in range(n_r_bins):
                df["{}_{}".format(run_dir,ir)] = bin_map[:,ir]
            valid_rundirs.append(run_dir)

        except OSError as e:
            print(run_dir," not found")
            continue
    
    file_name = "data/{}rad_rad_phi_{}.hdf5".format(mode,run_id)
    df.to_hdf(file_name,"vrad_vrad_phi",mode='w')
    with h5py.File(file_name,'r+') as f:
        f.create_dataset("valid_rundirs",data=np.array(valid_rundirs,dtype=bytes))


