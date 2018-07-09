import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

if __name__ == '__main__':
    # default parameters
    sfr_split = None
    plotDWind = True
    plotSFR = True
    plotAngle = True
    landscape = False
    
    skiprate = 1

#     run_id = "2014"
#     run_ids     = ["2014","2014","2015","2014","2015","2014"]
#     runs = ["q2redo","q2edd05redo","q2edd10_aniso1fixed","q2edd10redo","q2edd10_aniso3fixed","q2edd20redo"]
# #     edd = [.01,.05,.1,.1,.1,.2] # incorrect anyway
#     names = ["A","B","C1","C2","C3","D"]

#     run_ids     = ["1055","1055"]
#     runs = ["q2_SN","q2_SN_slow"]
#     names = ["ASS","AS"]
#     run_ids     = ["1058","1058"]
#     runs = ["restest0m02","restest0m04"]
#     names = ["0m02","0m04"]
#     outfile = "../figures/timeseries_summary_res_fast.pdf"

#     run_ids     = ["2014","2014","2015","2014","2015","2014","1055","1055","2018"]
#     runs = ["q2redo","q2edd05redo","q2edd10_aniso1fixed","q2edd10redo","q2edd10_aniso3fixed","q2edd20redo","q2_SN_slow","q2_SN","treetest"]
# #     edd = [.01,.05,.1,.1,.1,.2] # incorrect anyway
#     names = ["A","B","C1","C2","C3","D",r"A$_*$",r"A$_{**}$","X"]
#     sfr_split = [["A","B","D","X"],[r"A$_*$",r"A$_{**}$"],["C1","C2","C3"]]
#     outfile = "../figures/timeseries_summary_q2.pdf"

#     run_ids     = ["2014","2014","2015","2014","2015","2014","1055","1055"]
#     runs = ["q2redo","q2edd05redo","q2edd10_aniso1fixed","q2edd10redo","q2edd10_aniso3fixed","q2edd20redo","q2_SN_slow","q2_SN"]
#     names = ["A","B","C1","C2","C3","D",r"A$_*$",r"A$_{**}$"]
#     sfr_split = [["A","B","D","X"],[r"A$_*$",r"A$_{**}$"],["C1","C2","C3"]]
#     plotDWind = False
#     plotSFR = True
#     plotAngle = False
#     landscape = True
#     outfile = "../figures/sfr_summary_q2.pdf"

#     run_ids     = ["2014","2014","2014","2014","1055","1055"]
#     runs = ["q2redo","q2edd05redo","q2edd10redo","q2edd20redo","q2_SN_slow","q2_SN"]
#     names = [r"Lowest $L$",r"Low $L$",r"Medium $L$, Standard Anisotropy","Large $L$",r"Lowest $L$, low SNR",r"Lowest $L$, high SNR"]
#     sfr_split = None #[["A","B","D","X"],[r"A$_*$",r"A$_{**}$"],["C1","C2","C3"]]
#     plotDWind = False
#     plotSFR = True
#     plotAngle = False
#     landscape = True
#     outfile = "../figures/sfr_summary_q2.pdf"

#     run_ids     = ["2014","2014","2015","2014","2015","2014"]
#     runs = ["q2redo","q2edd05redo","q2edd10_aniso1fixed","q2edd10redo","q2edd10_aniso3fixed","q2edd20redo"]
#     names = [r"Lowest $L$",r"Low $L$",r"Medium $L$, Low Anisotropy",r"Medium $L$, Standard Anisotropy","Medium $L$, High Anisotropy","Large $L$"]
#     sfr_split = None
#     plotDWind = False
#     plotSFR = False
#     plotAngle = True
#     landscape = True
#     outfile = "../figures/wind_summary_aniso.pdf"

    

    names = ["02","04","06","08","1"]
#     run_ids = ["2019"]*len(names)
    run_ids = ["2020","2020","1060","1060","1060"]
    runs = ["restest0m"+x for x in names]
    outfile = "../figures/timeseries_summary_res_fast.pdf"
    
    if sfr_split is not None:
        nsfr_plots = len(sfr_split)
    else:
        nsfr_plots = 1
    nsp=0
    if plotDWind: nsp+=1
    if plotSFR:
        sf_sp_offset = nsp
        nsp+=nsfr_plots
    if plotAngle:
        angle_sp_offset = nsp
        nsp+=1
#     nsp = 2+nsfr_plots
    
#     nsp = 5
#     nsp = 3

    if landscape:
#         fig,sp = P.subplots(1,nsp,figsize=(12,3.),sharex=True)
        fig,sp = P.subplots(1,nsp,figsize=(3.,3.),sharex=True)
    else:
        fig,sp = P.subplots(nsp,1,figsize=(5,2.*nsp),sharex=True)
    if not isinstance(sp,np.ndarray):
        sp = np.array([sp,None])
    
    linewidth = 1.5

    for irun,(run_id,output_dir) in enumerate(zip(run_ids,runs)):
        print("Plotting: ",output_dir)
        
        if plotSFR:
            sf_data = np.loadtxt("data/sf"+run_id+output_dir+".dat")
            sf_sp = sp[sf_sp_offset]
            if sfr_split is not None:
                for isplit,sfr_plotlist in enumerate(sfr_split):
                    if names[irun] in sfr_plotlist:
                        sf_sp  = sp[isplit+sf_sp_offset]
            sf_sp.plot(sf_data[:,0],sf_data[:,2],'C'+str(irun),lw=linewidth,label=names[irun])
        
        wind_data = np.loadtxt("data/windangle_evolution"+run_id+output_dir+".dat")
        good_slice = (wind_data[:,4]>0.)
        if plotAngle:
            sp[angle_sp_offset].plot(wind_data[good_slice,1][::skiprate],wind_data[good_slice,2][::skiprate],lw=linewidth,label=names[irun])
        if plotDWind:
            sp[0].plot(wind_data[good_slice,1][::skiprate],wind_data[good_slice,4][::skiprate],label=names[irun],lw=linewidth)
#         sp[3].plot(wind_data[:,1],wind_data[:,4]/edd[irun])
    
    if plotSFR:
        for isp in range(sf_sp_offset,sf_sp_offset+nsfr_plots):
            sp[isp].set_ylabel(r"SFR (M$_\odot$/yr)")
            sp[isp].set_yscale('log')
            sp[isp].set_ylim([1.e-5,1.])
            sp[isp].set_yticks([1e-5,1e-4,1e-3,1e-2,1e-1])

    if plotAngle:
        sp[angle_sp_offset].set_ylabel(r"Wind angle ($^\circ$)")

    if plotDWind:
        sp[0].set_ylabel(r"$\dotM_w$ (M$_\odot$/yr)")
        sp[0].set_yscale('log')

#     sp[3].set_ylabel(r"$\dotM_w/f_\mathrm{edd}$ (M$_\odot$/yr)")
#     sp[3].set_yscale('log')
    
    for this_sp in sp[0:nsp]:
#         this_sp.set_xlim([0.,1.2])
        this_sp.set_xlim([0.,0.2])
    
    if plotDWind:
        sp[0].set_ylim([1.e-5,.2])
#     sp[1+nsfr_plots].set_ylim([0.,14.])
    if plotAngle:
        sp[angle_sp_offset].set_ylim([0.,16.])

    if landscape:
        for this_sp in sp[0:nsp]:
            this_sp.set_xlabel("$t$ (Myr)")
    else:
        sp[nsp-1].set_xlabel("$t$ (Myr)")
    if plotDWind:
        sp[0].legend(loc='best',fontsize='xx-small')
    else:
        for isp in range(0,nsfr_plots):
            sp[isp].legend(loc='best',fontsize='xx-small')

#     sp[0].legend(loc='best')
    
    P.tight_layout()
    P.savefig(outfile)
    P.close()