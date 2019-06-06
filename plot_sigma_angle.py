print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as P
import itertools

intrinsic_tau_transition = 40.
intrinsic_tau0 = 5.
intrinsic_index0 = 10.
intrinsic_index1 = 50.5

def tau_intrinsic(theta):
    if type(theta)==np.ndarray:
        below_slice = (theta<intrinsic_tau_transition)
        tau_out = np.zeros_like(theta)
        if np.sum(below_slice)>0:
            tau_out[below_slice]=intrinsic_tau0-theta[below_slice]/intrinsic_index0
        above_slice = ~below_slice
        if np.sum(above_slice)>0:
            tau_out[above_slice]=intrinsic_tau0-intrinsic_tau_transition/intrinsic_index0-(theta[above_slice]-intrinsic_tau_transition)/intrinsic_index1
        return tau_out
    elif type(theta)==float:
        if theta<intrinsic_tau_transition:
            return intrinsic_tau0-theta/intrinsic_index0
        else:
            return intrinsic_tau0-intrinsic_tau_transition/intrinsic_index0-(theta-intrinsic_tau_transition)/intrinsic_index1


# titlegroups = [[r"$f_{edd}=0.01$, $f_{a}=10^{2}$",r"$f_{edd}=0.05$, $f_{a}=10^{2}$",r"$f_{edd}=0.1$, $f_{a}=10^{2}$",r"$f_{edd}=0.2$, $f_{a}=10^{2}$"],[r"$f_{edd}=0.1$, $f_{a}=10^{1}$",r"$f_{edd}=0.1$, $f_{a}=10^{2}$",r"$f_{edd}=0.1$, $f_{a}=10^{3}$"]]
# titlegroups = [["A","B","C2","D"],["C1","C2","C3"],["A",r"A$_*$",r"A$_{**}$"],["A",r"A$_*$",r"A$_{**}$"]]
# colourgroups = [[0,1,3,5],[2,3,4],[0,6,7],[0,6,7]]
# irungroups = ["sigma_angle_many_edd_","sigma_angle_many_aniso_","sigma_angle_many_sf_","sigma_angle_many_sf_"]
# angle_maxes = [45.,45.,45.,90.]

# titlegroups = [["No SNR","Low SNR","High SNR"]]
# colourgroups = [[0,6,7]]
# irungroups = ["sigma_angle_many_sf_"]
# angle_maxes = [90.]

# for paper 1
# titlegroups = [["run_a0_e1","run_a1_e1","run_a2_e1","run_a3_e1"],["run_a2_e01","run_a2_e02","run_a2_e05","run_a2_e1","run_a2_e2"],["run_a2_e01_T1000","run_a2_e01_T300","run_a2_e01_T30","run_a2_e01"],["run_a2_e01","run_a2_e01_SN100","run_a2_e01_SN1000"]]
# titlegroups = [["a0_e1","a1_e1","a2_e1","a3_e1"],["a2_e01","a2_e02","a2_e05","a2_e1","a2_e2"],["a2_e01_T1000","a2_e01_T300","a2_e01_T30","a2_e01"],["a2_e01","a2_e01_SN100","a2_e01_SN1000"]]
# titlegroups = [["a0_e1","a1_e1","a2_e1","a3_e1"],["a2_e01","a2_e05","a2_e1","a2_e2"],["a2_e01_T1000","a2_e01_T300","a2_e01_T30","a2_e01"],["a2_e01","a2_e01_SN100","a2_e01_SN1000"]]
# colourgroups = [[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2]]
# irungroups = ["prodrun_aniso","prodrun_edd","prodrun_mintemp","prodrun_SN"]
# angle_maxes = [45.,45.,45.,45.]
# itimes = [20,100]
# scale = 1.


# for TORUS2018 presentation
# titlegroups = [["a0_e1","a1_e1","a2_e1","a3_e1"],["a2_e01","a2_e01_SN100","a2_e01_SN1000"]]
# colourgroups = [[0,1,2,3],[0,1,2]]
# irungroups = ["prodrun_aniso","prodrun_SN"]
# angle_maxes = [45.,45.,45.,45.]
# itimes = [100]
# scale = .5

# for Q=2 test
# titlegroups = [["SF_test_high_rho_floor_30_thinner_Q2_cut"]]
# colourgroups = [[0]]
# irungroups = ["prodrun_Q2_basis"]
# angle_maxes = [90.]
# itimes = [25,50,75,100]
# scale = .5
# transpose = False
# ncols = 4
# plot_tau = False

titlegroups = [  ["settled","vesc","rapid"],
                ["settled","vesc","rapid"]]
colourgroups = [[0,1,2],[0,1,2]]
irungroups = ["prodrun_equatorial","prodrun_polar"]
angle_maxes = [90.,90.]
itimes = [100]
scale = .5
transpose = True
ncols = 6
plot_tau = True
singleFile = "../figures/sigma_angle_longrun_weakflow_tau.pdf" 
# plot_tau = False
# singleFile = "../figures/sigma_angle_longrun_weakflow.pdf" 

# itimes = [100,200,500,1000]
# itimes = [1000]
# itimes = [200]

# titlegroups = [["blobrot","blobclose"]]
# colourgroups = [[0,1]]
# irungroups = ["prodrun_blob"]
# angle_maxes = [90.]
# 
# itimes = [100]


# time_between_dump = 1.e-7
time_between_dump = 1.e-5
conversion_to_Myr= 0.9778e9/1.e6

times = np.array(itimes)*time_between_dump*conversion_to_Myr

# singleFile = None # print out one figure per input
# singleFile = "../figures/sigma_angle_many_many_2014.pdf"
# singleFile = "../figures/sigma_angle_many_sf.pdf"
# singleFile = "../figures/sigma_angle_newruns.pdf"
# singleFile = "../figures/sigma_angle_blob.pdf"
# singleFile = "../figures/sigma_angle_newruns.pdf" # for dec 2018 paper
# singleFile = "../figures/sigma_angle_torus_2018.pdf" # for TORUS2018
# singleFile = "../figures/sigma_angle_Q2_basis.pdf" 

if plot_tau:
    col_offset = 3
else:
    col_offset = 0

if singleFile:
#     fig,subplorts = P.subplots(len(itimes),len(irungroups),sharey=True,figsize=(12,len(itimes)*3.))
#     fig,subplorts = P.subplots(len(itimes),len(irungroups),sharey=True,figsize=(12.*scale,6.*scale))
    ncol = len(itimes)
    nrow = len(irungroups)
    if transpose:
        ncol,nrow = nrow,ncol

    fig,subplorts = P.subplots(nrow,ncol,sharey=True,figsize=(6.*scale*ncol,6.*scale*nrow))        
    if isinstance(subplorts,np.ndarray):
        if subplorts.ndim==2:
            sp = subplorts
        else:
            sp = np.empty((nrow,ncol),dtype=object)
            if nrow==1:
                sp[0,:] = subplorts
            else:
                sp[:,0] = subplorts
    else:
        sp = np.empty((1,1),dtype=object)
        sp[0,0] = subplorts

for time_index,itime in enumerate(itimes):
    for irungroup,prefix in enumerate(irungroups):

        infile = "data/"+prefix+"%03d.dat"%itime
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
#         colors = P.rcParams['axes.color_cycle']
        colors = P.rcParams['axes.prop_cycle'].by_key()['color']
        dimcolors = [x+"60" for x in colors]
        
        if singleFile:
            if transpose:
                cur_sp = sp[time_index,irungroup]
            else:
                cur_sp = sp[irungroup,time_index]
        else:
            fig,cur_sp = P.subplots(1,1)
        

        cur_sp.set_xlim([0.,angle_maxes[irungroup]])

        if plot_tau:
            tau_interior = tau_intrinsic(degs)

        for icol,title in enumerate(titles):
        #    if icol in [0,1,3]:
            denses = data[:,icol*ncols+1+col_offset]
            lows = data[:,icol*ncols+2+col_offset]
            highs = data[:,icol*ncols+3+col_offset]
            #P.plot(degs,denses,color='black')
            # dim for intrinsic version
            if plot_tau:
                cur_sp.plot(degs,denses+tau_interior,color=dimcolors[colourgroups[irungroup][icol]],ls='--')
                cur_sp.errorbar(degs,denses+tau_interior,yerr=[denses-lows,highs-denses],color=dimcolors[colourgroups[irungroup][icol]],capsize=2)

            cur_sp.plot(degs,denses,color=colors[colourgroups[irungroup][icol]],label=title)
#             P.fill_between(degs,lows,highs,facecolor=colors[icol],alpha=.5)
            cur_sp.errorbar(degs,denses,yerr=[denses-lows,highs-denses],color=colors[colourgroups[irungroup][icol]],capsize=2)
#             print(title,time_index,icol,colourgroups[time_index][icol],colors[colourgroups[irungroup][icol]])

#         cur_sp.plot(degs,tau_intrinsic(degs),color='k',label="instrinsic")

        cur_sp.legend(fontsize='x-small')
#         cur_sp.set_ylim([1.e16,1.e27])
        if time_index==len(itimes)-1:
            cur_sp.set_xlabel(r"$\phi$ ($^\degree$)")
        #P.ylabel(r"$\Sigma$ ($M_\odot$/pc$^2$)")
        if plot_tau:
            cur_sp.set_ylim([1.e-5,5.])
            cur_sp.axhline(1.,c='k',ls='--')
        else:
            cur_sp.set_ylim([1.e15,1.e25])
            #hline for fun
            cur_sp.axhline(1.e22,c='k',ls='--')
        cur_sp.set_yscale('log')
        if irungroup==0:
            if plot_tau:
                cur_sp.set_ylabel(r"$N_\mathrm{H}$ (cm$^{-2}$)")
            else:
                cur_sp.set_ylabel(r"$\tau$")
        elif irungroup==len(irungroups)-1:
            cur_sp.yaxis.set_label_position("right")
            cur_sp.set_ylabel(r"t=%5.2f kyr"%(time*1.e3),size='x-large')

#         P.savefig(outp,dpi=200)
        
        
        if not singleFile:
            fig.savefig(outp)
            P.close()

if singleFile:
    fig.tight_layout()
    fig.savefig(singleFile)
    P.close()

