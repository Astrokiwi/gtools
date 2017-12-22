print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')

import pylab as P

#infile = "data/sigma_angle_sfnograv.dat"
#infile = "data/sigma_angle_sdev.out"
# infile = "data/aniso_sigma_angle.dat"
#infile = "data/aniso_sigma_angle_lowres.dat"
#titles = ["Grav","Grav+wSF","Grav+SF","NoGrav"]
#titles = ["SelfGrav","SelfGrav+SF","SelfGrav+SF+HighEdd","NoSelfGrav"]
# titles = ["SelfGrav","SelfGrav+SF","NoSelfGrav"]
#outp = "../figures/sigma_angle_sfnograv.png"
#outp = "../figures/sigma_angle_sfnograv_new.png"
# outp = "../figures/sigma_angle_sfnograv_aniso.png"
#outp = "../figures/sigma_angle_sfnograv_aniso_lowres.png"

# infile = "data/sigangle_small_big.dat"
# titles = ["small_big"]
# outp = "../figures/sigma_angle_small_big.png"

infile = "data/mondaytalk_aniso.dat"
titles = [r"$f_r/f_z=1$",r"$f_r/f_z=10$",r"$f_r/f_z=100$",r"$f_r/f_z=1000$"]
outp = "../figures/mondaytalk_aniso.png"

# infile = "data/mondaytalk_Edd.dat"
# titles = [r"$\epsilon=0.05$",r"$\epsilon=0.1$",r"$\epsilon=0.2$"]
# outp = "../figures/mondaytalk_Edd.png"

# infile = "data/mondaytalk_inclined.dat"
# titles = [r"$0^\circ$",r"$1^\circ$",r"$10^\circ$",r"$30^\circ$"]
# outp = "../figures/mondaytalk_inclined.png"


data = np.loadtxt(infile)

degs = data[:,0]/(2.*np.pi)*360

#colors=['brown','blue','yellow','purple']
colors = ['#0072b2', '#009e73', '#d55e00', '#56b4e9']


P.figure()
for icol,title in enumerate(titles):
#    if icol in [0,1,3]:
    denses = data[:,icol*3+1]
    lows = data[:,icol*3+2]
    highs = data[:,icol*3+3]
    #P.plot(degs,denses,color='black')
    P.plot(degs,denses,color=colors[icol],label=title)
    #P.fill_between(degs,lows,highs,facecolor=colors[icol],alpha=.5)

P.legend()
P.yscale('log')
P.xlabel(r"$\phi$ ($^\degree$)")
#P.ylabel(r"$\Sigma$ ($M_\odot$/pc$^2$)")
P.ylabel(r"$\Sigma$ (cm$^{-2}$)")
P.savefig(outp,dpi=200)
P.close()
