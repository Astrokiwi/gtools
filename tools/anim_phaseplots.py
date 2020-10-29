print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
# this_dir, this_filename = os.path.split(__file__)

#from joblib import Parallel, delayed
from multiprocessing import Pool
from functools import partial

from . import phaseplots
import re

from . import gizmo_tools

import glob

def sort_nicely( l ):
    """ Sort the given list in the way that humans expect.
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )


def g(run_id,output_dir,includedVal,rcut,m_bh,gridMode,uniqueMode,snap_str):
    return phaseplots.savephaseplots(run_id, output_dir, snap_str, includedVal, rcut=rcut, m_bh=m_bh, gridMode=gridMode, uniqueMode=uniqueMode)


if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]

    try:
        nprocs = int(sys.argv[3])
    except IndexError:
        nprocs = 64 # default
    

    print("Running")

#     gridMode = False
#     uniqueMode = False
    gridMode = True
    uniqueMode = False
    
#     append_name = "_high_hcn"
    append_name = ""

    #includedVals = ["dt_p","nH_p","TK_p","rad2d_p","z_p","vrad","dHeat","dt_heat"]
#     includedVals = ["rad_p","nH_p","arad_p","radrad_p"]
#     includedVals = ["rad2d_p","radrad_p"]
#     includedVals = ["co1","co2"]
#     includedVals = ["hcn2","nH_p","TK_p","tau","flux_p"]
#     includedVals = ["nH_p","TK_p","tau","flux_p"]
#     includedVals = ["nH_p","tau"]
    includedVals = ["nH_p","TK_p"]
#     includedVals = ["vrad","TK_p"]
#     includedVals = ["nH_p","arad_p"]
#     includedVals = ["arad_p","nH_p","TK_p","flux_p","rad_p","prad","pratio"]
#     includedVals = ["nH_p","TK_p","vrad","rad_p","mJ_p"]
#     includedVals = ["nH_p","age"]
#     includedVals = ["ylos","vlos"]
#     includedVals = ["xlos","vlos","vrad","vcirc"]
#     includedVals = ["ylos","vlos","vrad","vcirc"]
#     includedVals = ["xlos","ylos"]
#     includedVals = ["vcirc","vrad"]#,"tau"]
#     includedVals = ["phi","vrad"]
#     includedVals = ["phi","vcirc"]
#     includedVals = ["vrad","TK_p",""]
#     includedVals = ["TK_p","vel"]
#     includedVals = ["tau","pratio"]
#     includedVals = ["p_p","prad"]
#     includedVals = ["rad_p","vrad"]
#     includedVals = ["hcn2","nH_p"]
#     includedVals = ["hcn1","hcn2","co1","co2"]
#     includedVals = ["rad_p","radrad_p","nH_p"]
#     includedVals = ["dustTemp","TK_p","nH_p"]
    #includedVals = ["rad2d_p","hz_rat"]
    #includedVals = ["mJ_p","nH_p","TK_p"]
#     includedVals = ["arad_p","radrad_p"]
#     includedVals = ["vrad","z_p"]
    #includedVals = ["flux_p","radrad_p"]
#     includedVals = ["nH_p","TK_p"]
#     includedVals = ["cs_p","vrad"]
#     includedVals = ["rad2d_p","agn_heat_p","radrad_p"]
#     includedVals = ["TK_p","opac"]
    #includedVals = ["p_p","nH_p"]
    #includedVals = ["rad2d_p","nH_p"]
#     includedVals = ["dt_p","rad_p","nH_p","TK_p"]
    #includedVals = ["depth_p","radrad_p"]
#     includedVals = ["dt_p","TK_p"]
#     includedVals = ["TK_p","agn_heat_p"]
#     includedVals = ["vel","TK_p"]
#     includedVals = ["rad_p","nH_p","TK_p","opac","arad_p","radrad_p"]
#     includedVals = ["rad_p","arad_p"]
    #includedVals = ["prat"]
    #includedVals = ["nH_p","p_p"]
    #includedVals = ["h_p","nH_p"]
    #includedVals = ["mJ_p","prat"]
#     includedVals = ["vcirc","vrad"]
    
    
#     rcut = 80.
    rcut = 1.e9
#     rcut = 2.

    m_bh = 1.e6
#     m_bh = .05e6

#     snapi = 0

    snapi = 0
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir,False)
#     snapf = 22
#     print("HACKING SNAPF to 22")
    
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
    varOneRange = [0]
    if gridMode:
        varOneRange = range(l)
#     if uniqueMode:
#         varOneRange = range(0,l,2)
    for i in varOneRange: # do first vs everything
        if uniqueMode and i%2!=0:
            continue
        varTwoRange = range(i+1,l)
        if uniqueMode:
            varTwoRange = [i+1]
        for j in varTwoRange:
            anim_names.append(includedVals[i]+includedVals[j])
            n_anim+=1
    
    pool = Pool(processes=nprocs)
    h = partial(g,run_id,output_dir,includedVals,rcut,m_bh,gridMode,uniqueMode)
#     animstrs = pool.starmap(phaseplots.savephaseplots,zip(it.repeat(run_id),it.repeat(output_dir),snap_strs,it.repeat(includedVals),it.repeat(rcut)))
    animstrs = pool.map(h,snap_strs)
    pool.close()
    animstrs = np.vstack(animstrs)

    for i_anim in range(n_anim):
        #anim_file_list = "data/phaseplot_animlist"+run_id+"_"+output_dir+"_"+anim_names[i_anim]+".txt"
        #np.savetxt(anim_file_list,animstrs[:,i_anim],fmt='%s')
        first_filename = animstrs[0,i_anim]
        anim_catch = first_filename[:-7]+"%03d.png"
        cmd = "ffmpeg -y -r 24 -i "+anim_catch+" -c:v mpeg4 -q:v 1 "+movieDir+"/phaseplots"+run_id+"_"+output_dir+"_"+anim_names[i_anim]+append_name+".mp4"
        # remove all old images
        n_last = animstrs.shape[0]
        for imageFile in glob.glob(first_filename[:-7]+"*.png"):
            if int(imageFile[-7:-4])>=n_last:
                os.remove(imageFile)
        
        #cmd = "ffmpeg -y -r 5 -i "+anim_file_list+" -c:v mpeg4 -q:v 1 /export/1/djw/movies/phaseplots"+run_id+"_"+output_dir+"_"+anim_names[i_anim]+".mp4"
        os.system(cmd)
# 
#     for i in range(snapf):
#         x = phaseplots.savephaseplots(run_id,output_dir,snap_strs[i],includedVals)
#         animstrs[i] = x
#         print(animstrs)
    
    
    
    
