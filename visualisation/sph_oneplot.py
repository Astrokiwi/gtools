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

run_id = sys.argv[1]
output_dir = sys.argv[2]

snapx = int(sys.argv[3])


gizmoDir = gizmo_tools.getGizmoDir()
movieDir = gizmo_tools.getMovieDir()
fullDir = gizmoDir+"/"+run_id+"/"+output_dir


infile = fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5"
outfile = "../pics/sphoneplot"+run_id+output_dir+"%03d.png"%snapx
#outfile = "../pics/georgia_sphoneplot"+run_id+output_dir+"%03d.png"%snapx
#outfile = "../pics/pretty_sphoneplot"+run_id+output_dir+"%03d.png"%snapx
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',L=512,flat=True,ring=True,nrows=1)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['vels','dens'],L=400,scale=20.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['vels'],L=600,scale=20.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens'],L=800,scale=40.,cols=1)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['dens'],views=['side'],L=400,scale=2.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['view'],L=400,scale=2.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='prism',flat=False,ring=False,vorinoi=True,plot=['rand'],L=1200,scale=30.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',views=['face'],flat=True,ring=True,plot=['view'],L=800,scale=3.,rot=[0,80./360.*2.*np.pi])


# outfile = "../pics/sphviewplot"+run_id+output_dir+"%03d.png"%snapx
# sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',views=['face'],flat=True,ring=True,plot=['view'],L=800,scale=1.,rot=[0,80./360.*2.*np.pi])
#outfile = "../pics/sphmanyplot"+run_id+output_dir+"%03d.png"%snapx
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',views=['face'],flat=True,ring=True,plot=['tdust'],L=800,scale=1.,rot=[0,80./360.*2.*np.pi])
# outfile = "../pics/sphmaxplot"+run_id+output_dir+"%03d.png"%snapx
# sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',views=['face'],flat=True,ring=True,plot=['nH','temp','opac','tdust'],L=800,scale=1.,rot=[0,80./360.*2.*np.pi],maxmode=True)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',views=['face'],flat=True,ring=True,plot=['nH'],L=800,scale=1.,rot=[0,80./360.*2.*np.pi],maxmode=True)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',views=['face'],flat=True,ring=True,plot=['nH'],L=800,scale=2.)
sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',views=['face'],flat=True,ring=True,plot=['temp'],L=800,scale=1.e4)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['dens','temp'],L=400,scale=2.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens','temp'],L=800,scale=5.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens'],L=600)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens'],L=400,planenorm=True,scale=10.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['dens','col'],L=400,scale=10.)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='Spectral',flat=True,ring=True,plot=['vel_2d','vel_x','vel_y','vel_z','vel_r'],L=600)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='Spectral',flat=False,ring=False,plot=['vel_a'],L=400,scale=2.)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['dens'],L=800,scale=8.)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['arad','col'],L=800,scale=1.5)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['dens'],L=800,scale=1.5)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['tau'],L=800,scale=3)


#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['vel_2d','vel_z','vels'],L=400,scale=1.)
#,views=['side']
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['dens'],L=800,scale=4.,views=['face'])
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['temp'],L=1024,scale=10.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['temp'],L=400,subsample=10.,pixsize=1)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='Greys',flat=True,ring=True,plot=['emit','temp','dens'],L=400)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=False,ring=False,plot=['view'],scale=3.,L=600,rot=[0.,np.pi/4.])

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='coolwarm',flat=False,ring=False,plot=['vlos'],scale=40.,L=600)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=False,plot=['dens'],L=600,cols=1,scale=15.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=False,plot=['dens'],L=600,cols=1,scale=15.,visibleAxes=False)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens'],L=10000,views=['face','side'],scale=10.,visibleAxes=False)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['dt'],L=200,views=['face','side'],scale=15.,pixsize=2)


#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dust','dg','dens'])
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens','temp','vels'],L=200,scale=10.)
#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['dens','temp','vels'],views=['side'],L=200,scale=10.)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['dens'],L=100,cols=1,visibleAxes=False)

#sph_frame.makesph_trhoz_frame(infile,outfile,cmap='plasma',flat=True,ring=True,plot=['col'])

# outfile = "../pics/georgia_sphvelplot"+run_id+output_dir+"%03d.png"%snapx
# sph_frame.makesph_trhoz_frame(infile,outfile,cmap='viridis',flat=True,ring=True,plot=['vels'],L=400,scale=10.)
# outfile = "../pics/georgia_sphzonesplot"+run_id+output_dir+"%03d.png"%snapx
# sph_frame.makesph_trhoz_frame(infile,outfile,cmap='RGBsteps',flat=True,ring=True,plot=['nH'],L=400,scale=10.)

