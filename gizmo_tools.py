import socket
import os
import numpy as np
import re

molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16
msun_g = 1.989e33

def sort_nicely( l ):
    """ Sort the given list in the way that humans expect.
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )

# def getGizmoDir():
#     sname = socket.gethostname()
# 
#     if ( sname=="trillian" ):
#         return "/export/1/djw/gizmos"
#     elif ( sname=="srv01921" ):
#         return "/srv/djw1g16/gizmos"
#     
#     raise Exception("Unknown server; add server and directory to gizmodatadir.py")


def getGizmoDir(irun):
    sname = socket.gethostname()

    if ( sname=="trillian" ):
        paths = ["/export/1/djw/gizmos","/export/2/djw/gizmos"]
        for path in paths:
            if os.path.isdir("{}/{}".format(path,irun)):
                return path
#         return "/export/1/djw/gizmos"
        raise Exception("Directory not found on {}".format(sname))
    elif ( sname=="srv01921" ):
        return "/srv/djw1g16/gizmos"
    
    raise Exception("Unknown server; add server and directory to gizmodatadir.py")

def getMovieDir():
    sname = socket.gethostname()

    if ( sname=="trillian" ):
        return "/export/1/djw/movies"
    elif ( sname=="srv01921" ):
        return "/srv/djw1g16/movies"
    
    raise Exception("Unknown server; add server and directory to gizmodatadir.py")

def lastConsecutiveSnapshot(run_id,output_dir):
    gizmoDir = getGizmoDir(run_id)
    movieDir = getMovieDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir

    fnames = os.listdir(fullDir)
    sort_nicely(fnames)
    fnames = np.array(fnames)
    #fnames.sort()
    snapshotfilebools = np.array([x.startswith("snapshot") for x in fnames])
    snapshotfiles = fnames[snapshotfilebools]
    
#     print(snapshotfiles)
    
    snapf = 0

    ctime = os.path.getmtime(fullDir+"/snapshot_000.hdf5")
    for fname in snapshotfiles[1:]:
        new_snapf = int(fname[9:len(fname)-5])
        new_ctime = os.path.getmtime(fullDir+"/"+fname)
#         print(ctime,new_ctime)
        if new_ctime>ctime and snapf==new_snapf-1 :
            ctime = new_ctime
            snapf = new_snapf
    
    return snapf