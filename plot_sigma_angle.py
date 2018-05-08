print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as P

# titlegroups = [[r"$f_{edd}=0.01$, $f_{a}=10^{2}$",r"$f_{edd}=0.05$, $f_{a}=10^{2}$",r"$f_{edd}=0.1$, $f_{a}=10^{2}$",r"$f_{edd}=0.2$, $f_{a}=10^{2}$"],[r"$f_{edd}=0.1$, $f_{a}=10^{1}$",r"$f_{edd}=0.1$, $f_{a}=10^{2}$",r"$f_{edd}=0.1$, $f_{a}=10^{3}$"]]
titlegroups = [["A","B","C2","D"],["C1","C2","C3"]]
colourgroups = [[0,1,3,5],[2,3,4]]
irungroups = ["sigma_angle_many_edd_","sigma_angle_many_aniso_"]

# itimes = [100,200,500,1000]
# itimes = [1000]
itimes = [200]

time_between_dump = 1.e-6
conversion_to_Myr= 0.9778e9/1.e6

times = np.array(itimes)*time_between_dump*conversion_to_Myr

# singleFile = None # print out one figure per input
singleFile = "../figures/sigma_angle_many_many_2014.pdf"

if singleFile:
    fig,subplorts = P.subplots(len(itimes),len(irungroups),sharex=True,sharey=True,figsize=(12,len(itimes)*4.))
    if subplorts.ndim==2:
        sp = subplorts
    else:
        sp = np.empty((1,subplorts.shape[0]),dtype=object)
        sp[0,:] = subplorts

for time_index,itime in enumerate(itimes):
    for irungroup,prefix in enumerate(irungroups):

        infile = "data/"+prefix+str(itime)+".dat"
        titles = titlegroups[irungroup]
        time = times[time_index]
#         outp = "../figures/"+prefix+str(itime)+".png"
        if not singleFile:
            outp = "../figures/"+prefix+str(itime)+".pdf"

        print(infile,time)

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

        # infile = "data/sigma_angle_many_edd_500.dat"
        # titles = [r"$f_{edd}=0.01$",r"$f_{edd}=0.05$",r"$f_{edd}=0.1$",r"$f_{edd}=0.2$"]
        # outp = "../figures/sigma_angle_many_edd_500.png"

#         infile = "data/sigma_angle_many_aniso_500.dat"
#         titles = [r"$f_{a}=10^{1}$",r"$f_{a}=10^{2}$",r"$f_{a}=10^{2}$"]
#         outp = "../figures/sigma_angle_many_aniso_500.png"

        # infile = "data/mondaytalk_Edd.dat"
        # titles = [r"$\epsilon=0.05$",r"$\epsilon=0.1$",r"$\epsilon=0.2$"]
        # outp = "../figures/mondaytalk_Edd.png"

        # infile = "data/mondaytalk_inclined.dat"
        # titles = [r"$0^\circ$",r"$1^\circ$",r"$10^\circ$",r"$30^\circ$"]
        # outp = "../figures/mondaytalk_inclined.png"

        # infile = "data/mondaytalk_aniso.dat"
        # titles = [r"$f_r/f_z=1$",r"$f_r/f_z=10$",r"$f_r/f_z=100$",r"$f_r/f_z=1000$"]
        # outp = "../figures/mondaytalk_aniso.png"


        data = np.loadtxt(infile)

        degs = data[:,0]/(2.*np.pi)*360

        #colors=['brown','blue','yellow','purple']
#         colors = ['#0072b2', '#009e73', '#d55e00', '#56b4e9']
        colors = P.rcParams['axes.color_cycle']


        
        if singleFile:
            cur_sp = sp[time_index,irungroup]
        else:
            fig,cur_sp = P.subplots(1,1)
        

        cur_sp.set_xlim([0.,45.])

        for icol,title in enumerate(titles):
        #    if icol in [0,1,3]:
            denses = data[:,icol*4+1]
            #P.plot(degs,denses,color='black')
            cur_sp.plot(degs,denses,color=colors[colourgroups[irungroup][icol]],label=title)
            lows = data[:,icol*4+2]
            highs = data[:,icol*4+3]
#             P.fill_between(degs,lows,highs,facecolor=colors[icol],alpha=.5)
            cur_sp.errorbar(degs,denses,yerr=[denses-lows,highs-denses],color=colors[colourgroups[irungroup][icol]],capsize=2)
#             print(title,time_index,icol,colourgroups[time_index][icol],colors[colourgroups[irungroup][icol]])

        cur_sp.legend(fontsize='x-small')
        cur_sp.set_yscale('log')
        cur_sp.set_ylim([1.e11,1.e27])
        if time_index==len(itimes)-1:
            cur_sp.set_xlabel(r"$\phi$ ($^\degree$)")
        #P.ylabel(r"$\Sigma$ ($M_\odot$/pc$^2$)")
        if irungroup==0:
            cur_sp.set_ylabel(r"$N_\mathrm{H}$ (cm$^{-2}$)")
        if irungroup==len(irungroups)-1:
            cur_sp.yaxis.set_label_position("right")
            cur_sp.set_ylabel(r"t=%5.3f Myr"%time)

#         P.savefig(outp,dpi=200)
        
        if not singleFile:
            fig.savefig(outp)
            P.close()

if singleFile:
    fig.tight_layout()
    fig.savefig(singleFile)
    P.close()

