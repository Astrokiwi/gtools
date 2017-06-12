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

run_id = sys.argv[1]
output_dir = sys.argv[2]

snapx = int(sys.argv[3])

fullDir = "/export/1/djw/gizmos/"+run_id+"/"+output_dir

infile = fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5"
outfile = "../pics/sphoneplot"+run_id+output_dir+"%03d.png"%snapx
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',L=512,flat=True,ring=True,nrows=1)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['vels','dens'],L=400,subsample=10,pixsize=1)
sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens'],L=400)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['tdust','temp'],L=400,scale=10.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens','temp'],L=800,scale=5.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens','temp'],L=400)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens'],L=400,planenorm=True,scale=10.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['dens','col'],L=400,scale=10.)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['temp','dens'],L=400,scale=20.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['temp'],L=1024,scale=10.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['temp'],L=400,subsample=10.,pixsize=1)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='Greys',flat=True,ring=True,plot=['emit','temp','dens'],L=400)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['view'],L=400)
