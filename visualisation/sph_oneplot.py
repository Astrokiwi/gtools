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
#outfile = "../pics/georgia_sphoneplot"+run_id+output_dir+"%03d.png"%snapx
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',L=512,flat=True,ring=True,nrows=1)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['vels','dens'],L=400,scale=20.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['vels'],L=600,scale=20.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens'],L=800,scale=40.,cols=1)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['tdust','temp'],L=400,scale=10.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens','temp'],L=800,scale=5.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens','temp'],L=400)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens'],L=400,planenorm=True,scale=10.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['dens','col'],L=400,scale=10.)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['temp','dens'],L=400,scale=20.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['temp'],L=1024,scale=10.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['temp'],L=400,subsample=10.,pixsize=1)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='Greys',flat=True,ring=True,plot=['emit','temp','dens'],L=400)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['view'],scale=40.,L=600,cols=1,rot=[0.,np.pi/2.])
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='coolwarm',flat=False,ring=False,plot=['vlos'],scale=40.,L=600)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=False,plot=['dens'],L=600,cols=1,scale=15.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=False,plot=['dens'],L=600,cols=1,scale=15.,visibleAxes=False)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dust','dg','dens'])
sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens','temp','vels'],L=200,scale=10.)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['dens'],L=100,cols=1,visibleAxes=False)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['col'])

# outfile = "../pics/georgia_sphvelplot"+run_id+output_dir+"%03d.png"%snapx
# sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['vels'],L=400,scale=10.)
# outfile = "../pics/georgia_sphzonesplot"+run_id+output_dir+"%03d.png"%snapx
# sph_frame.makesph_trhoz_frame(infile,outfile,cmap='RGBsteps',flat=True,ring=True,plot=['nH'],L=400,scale=10.)

