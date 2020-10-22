print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
#from joblib import Parallel, delayed

#from sys import path
#path.append("../src/")
import sph_frame

print("Running")

#run_id = sys.argv[1]
#output_dir = sys.argv[2]

#snapx = int(sys.argv[3])

run_ids =       ["1016",        "1016",     "1016",     "1019"]
output_dirs =   ["fataradtest", "sftest",   "sfmild",   "nograv"]
snapx = 56

for run_id,output_dir in zip(run_ids,output_dirs):


    fullDir = "/export/1/djw/gizmos/"+run_id+"/"+output_dir
    infile = fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5"

    outfile = "../pics/nHdiscrete"+run_id+output_dir+"%03d.png"%snapx
    sph_frame.makesph_trhoz_frame(infile,outfile,cmap='RGBsteps',flat=True,ring=True,plot=['nH'],L=400,scale=20.)
    outfile = "../pics/nHsmooth"+run_id+output_dir+"%03d.png"%snapx
    sph_frame.makesph_trhoz_frame(infile,outfile,flat=True,ring=True,plot=['nH'],L=400,scale=20.)

