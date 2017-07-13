print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')

import pylab as P

#infile = "data/sigma_angle_sfnograv.dat"
#infile = "data/sigma_angle_sdev.out"
infile = "data/aniso_sigma_angle.dat"
#infile = "data/aniso_sigma_angle_lowres.dat"
#titles = ["Grav","Grav+wSF","Grav+SF","NoGrav"]
#titles = ["SelfGrav","SelfGrav+SF","SelfGrav+SF+HighEdd","NoSelfGrav"]
titles = ["SelfGrav","SelfGrav+SF","NoSelfGrav"]
#outp = "../figures/sigma_angle_sfnograv.png"
#outp = "../figures/sigma_angle_sfnograv_new.png"
outp = "../figures/sigma_angle_sfnograv_aniso.png"
#outp = "../figures/sigma_angle_sfnograv_aniso_lowres.png"

data = np.loadtxt(infile)

degs = data[:,0]/(2.*np.pi)*360

colors=['brown','blue','yellow'] #,'purple'

P.figure()
for icol,title in enumerate(titles):
#    if icol in [0,1,3]:
    denses = data[:,icol*3+1]
    lows = data[:,icol*3+2]
    highs = data[:,icol*3+3]
    P.plot(degs,denses,color='black')
    P.fill_between(degs,lows,highs,facecolor=colors[icol],alpha=.5,label=title)

P.legend()
P.yscale('log')
P.xlabel(r"$\phi$ ($^\degree$)")
#P.ylabel(r"$\Sigma$ ($M_\odot$/pc$^2$)")
P.ylabel(r"$\Sigma$ (cm$^{-2}$)")
P.savefig(outp,dpi=200)
P.close()
