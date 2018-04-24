print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
#from joblib import Parallel, delayed

#from sys import path
#path.append("../src/")
import sph_frame

from sys import path
path.append("../")
import gizmo_tools

print("Running")

# output_dir = "q2redo"
runs = ["q2edd05redo","q2edd10_aniso1redo","q2edd10_aniso3redo","q2edd10redo","q2edd20redo","q2redo"]
run_id = "2014"

# snapx = [0,100,200,500,1000]
snapx = 200

scales = [1.,4.,20.]
nscales = len(scales)

gizmoDir = gizmo_tools.getGizmoDir()
movieDir = gizmo_tools.getMovieDir()

for output_dir in runs:
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    infile = fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5"

    outfile = "../../figures/nHzoom"+run_id+output_dir+".png"
    sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['nH']*nscales,L=400,scale=scales)

    outfile = "../../figures/Tzoom"+run_id+output_dir+".png"
    sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['temp']*nscales,L=400,scale=scales)
