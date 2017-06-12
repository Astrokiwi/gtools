print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os

import re

import h5py



def sort_nicely( l ):
    """ Sort the given list in the way that humans expect.
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )

if __name__ == '__main__':
    print("Running")

    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    
    if ( len(sys.argv)>3 ):
        lossrad = float(sys.argv[3])
    else:
        lossrad = 20.

    snapi = 0
    fullDir = "/export/1/djw/gizmos/"+run_id+"/"+output_dir

    fnames = os.listdir(fullDir)
    sort_nicely(fnames)
    fnames = np.array(fnames)
    #fnames.sort()
    snapshotfilebools = np.array([x.startswith("snapshot") for x in fnames])
    snapshotfiles = fnames[snapshotfilebools]

    snapf = 0
    ctime = os.path.getmtime(fullDir+"/snapshot_000.hdf5")
    for fname in snapshotfiles[1:]:
        new_snapf = int(fname[9:len(fname)-5])
        new_ctime = os.path.getmtime(fullDir+"/"+fname)
        if ( new_ctime>ctime ) :
            ctime = new_ctime
            snapf = new_snapf

    os.system("rm ../pics/sphplot"+run_id+output_dir+"???.png")


    print("nfiles:",snapf-snapi+1)
    
    infiles = [fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5" for snapx in range(snapi,snapf+1)]

    dout = []
    for snapx,infile in enumerate(infiles):
        f = h5py.File(infile,"r")
    
        header = f["/Header"]
        time = header.attrs.get("Time")
        time*= 0.9778e9 # to yr
        time/=1.e6 # to Myr

        xyz_p = np.array(f["/PartType0/Coordinates"]) # kpc
        rad_p = np.sqrt(np.sum(xyz_p**2,axis=1))*1.e3 # to pc

        m_p = np.array(f["/PartType0/Masses"]) # 10^10 msun
        m_p*=1e+10 # 10^10 solar masses to solar masses
        
        inside = (rad_p<lossrad)
        m_inside = np.sum(m_p[inside])
        m_outside = np.sum(m_p[~inside])
        
        dout.append([time,m_inside,m_outside])
        print([time,m_inside,m_outside])
    
    np.savetxt("data/mrate_{}_{}.dat".format(run_id,output_dir),dout)
    
    
    
    
    
    
    