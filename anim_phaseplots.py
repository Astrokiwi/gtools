print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
#from joblib import Parallel, delayed
from multiprocessing import Pool
from functools import partial

import phaseplots
import re
import itertools as it

import gizmo_tools

import glob

def sort_nicely( l ):
    """ Sort the given list in the way that humans expect.
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )


def g(run_id,output_dir,includedVal,rcut,m_bh,gridMode,snap_str):
    return phaseplots.savephaseplots(run_id,output_dir,snap_str,includedVal,rcut=rcut,m_bh=m_bh,gridMode=gridMode)


if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]

    print("Running")

    gridMode = False

    #includedVals = ["dt_p","nH_p","TK_p","rad2d_p","z_p","vrad","dHeat","dt_heat"]
#     includedVals = ["rad_p","nH_p","arad_p","radrad_p"]
#     includedVals = ["rad2d_p","radrad_p"]
#     includedVals = ["co1","co2"]
    includedVals = ["hcn1","nH_p","TK_p","tau","flux_p"]
#     includedVals = ["rad_p","radrad_p","nH_p"]
#     includedVals = ["dustTemp","TK_p","vrad"]
    #includedVals = ["rad2d_p","hz_rat"]
    #includedVals = ["mJ_p","nH_p","TK_p"]
#     includedVals = ["arad_p","radrad_p"]
#     includedVals = ["vrad","z_p"]
    #includedVals = ["flux_p","radrad_p"]
#     includedVals = ["nH_p","TK_p"]
#     includedVals = ["rad2d_p","agn_heat_p","radrad_p"]
#     includedVals = ["TK_p","opac"]
    #includedVals = ["p_p","nH_p"]
    #includedVals = ["rad2d_p","nH_p"]
#     includedVals = ["dt_p","rad_p","nH_p","TK_p"]
    #includedVals = ["depth_p","radrad_p"]
#     includedVals = ["dt_p","TK_p"]
#     includedVals = ["TK_p","agn_heat_p"]
#     includedVals = ["rad_p","vel"]
#     includedVals = ["rad_p","nH_p","TK_p","opac","arad_p","radrad_p"]
#     includedVals = ["rad_p","arad_p"]
    #includedVals = ["prat"]
    #includedVals = ["nH_p","p_p"]
    #includedVals = ["h_p","nH_p"]
    #includedVals = ["mJ_p","prat"]
#     includedVals = ["vcirc","vrad"]
    
    
#     rcut = 80.
    rcut = 1.e9

    m_bh = 1.e6
#     m_bh = .05e6

#     snapi = 0

    snapi = 0
#     snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir)
    snapf = 100
    
    movieDir = gizmo_tools.getMovieDir()

#     fnames = os.listdir(fullDir)
#     sort_nicely(fnames)
#     fnames = np.array(fnames)
#     snapshotfilebools = np.array([x.startswith("snapshot") for x in fnames])
#     snapshotfiles = fnames[snapshotfilebools]
# 
#     snapf = 0
#     ctime = os.path.getmtime(fullDir+"/snapshot_000.hdf5")
#     for fname in snapshotfiles[1:]:
#         new_snapf = int(fname[9:len(fname)-5])
#         new_ctime = os.path.getmtime(fullDir+"/"+fname)
#         if ( new_ctime>ctime ) :
#             ctime = new_ctime
#             snapf = new_snapf

    
    
    l = len(includedVals)
#     n_anim = (l*(l-1))//2
    
    
    animstrs = np.empty((snapf),dtype=object)
    snap_strs = ["%03d" % i for i in range(snapf+1)]

    n_anim = 0
    anim_names = []
#     for i in range(l): # do full grid
    varOneRange = range(l) if gridMode else [0]
    for i in varOneRange: # do first vs everything
        for j in range(i+1,l):
            anim_names.append(includedVals[i]+includedVals[j])
            n_anim+=1
    
    pool = Pool(processes=80)
    h = partial(g,run_id,output_dir,includedVals,rcut,m_bh,gridMode)
#     animstrs = pool.starmap(phaseplots.savephaseplots,zip(it.repeat(run_id),it.repeat(output_dir),snap_strs,it.repeat(includedVals),it.repeat(rcut)))
    animstrs = pool.map(h,snap_strs)
    pool.close()
    animstrs = np.vstack(animstrs)

    for i_anim in range(n_anim):
        #anim_file_list = "data/phaseplot_animlist"+run_id+"_"+output_dir+"_"+anim_names[i_anim]+".txt"
        #np.savetxt(anim_file_list,animstrs[:,i_anim],fmt='%s')
        first_filename = animstrs[0,i_anim]
        anim_catch = first_filename[:-7]+"%03d.png"
        # remove all old images
        n_last = animstrs.shape[0]
        for imageFile in glob.glob(first_filename[:-7]+"*.png"):
            if int(imageFile[-7:-4])>=n_last:
                os.remove(imageFile)
        
        cmd = "ffmpeg -y -r 24 -i "+anim_catch+" -c:v mpeg4 -q:v 1 "+movieDir+"/phaseplots"+run_id+"_"+output_dir+"_"+anim_names[i_anim]+".mp4"
        #cmd = "ffmpeg -y -r 5 -i "+anim_file_list+" -c:v mpeg4 -q:v 1 /export/1/djw/movies/phaseplots"+run_id+"_"+output_dir+"_"+anim_names[i_anim]+".mp4"
        os.system(cmd)
# 
#     for i in range(snapf):
#         x = phaseplots.savephaseplots(run_id,output_dir,snap_strs[i],includedVals)
#         animstrs[i] = x
#         print(animstrs)
    
    
    
    
