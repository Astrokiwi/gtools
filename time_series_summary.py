import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

if __name__ == '__main__':
    sfr_split = None

#     run_id = "2014"
#     run_ids     = ["2014","2014","2015","2014","2015","2014","1055","1055","2018"]
#     runs = ["q2redo","q2edd05redo","q2edd10_aniso1fixed","q2edd10redo","q2edd10_aniso3fixed","q2edd20redo","q2_SN_slow","q2_SN","treetest"]
# #     edd = [.01,.05,.1,.1,.1,.2] # incorrect anyway
#     names = ["A","B","C1","C2","C3","D",r"A$_*$",r"A$_{**}$","X"]
#     sfr_split = [["A","B","D","X"],[r"A$_*$",r"A$_{**}$"],["C1","C2","C3"]]
#     outfile = "../figures/timeseries_summary_q2.pdf"

    names = ["02","04","06","08","1"]
    run_ids = ["2019"]*len(names)
    runs = ["restest0m"+x for x in names]
    outfile = "../figures/timeseries_summary_res.pdf"
    
    if sfr_split is not None:
        nsfr_plots = len(sfr_split)
    else:
        nsfr_plots = 1
    nsp = 2+nsfr_plots
    
#     nsp = 5
#     nsp = 3

    fig,sp = P.subplots(nsp,1,figsize=(5,8.),sharex=True)
    
    skiprate = 1
    
    linewidth = 1.5

    for irun,(run_id,output_dir) in enumerate(zip(run_ids,runs)):
        print("Plotting: ",output_dir)
        
        sf_data = np.loadtxt("data/sf"+run_id+output_dir+".dat")
        sf_sp = sp[1]
        if sfr_split is not None:
            for isplit,sfr_plotlist in enumerate(sfr_split):
                if names[irun] in sfr_plotlist:
                    sf_sp  = sp[isplit+1]
        sf_sp.plot(sf_data[:,0][::skiprate],sf_data[:,2][::skiprate],'C'+str(irun),lw=linewidth)
        
        wind_data = np.loadtxt("data/windangle_evolution"+run_id+output_dir+".dat")
        good_slice = (wind_data[:,4]>0.)
        sp[1+nsfr_plots].plot(wind_data[good_slice,1][::skiprate],wind_data[good_slice,2][::skiprate],lw=linewidth)
        sp[0].plot(wind_data[good_slice,1][::skiprate],wind_data[good_slice,4][::skiprate],label=names[irun],lw=linewidth)
#         sp[3].plot(wind_data[:,1],wind_data[:,4]/edd[irun])
    
    for isp in range(1,1+nsfr_plots):
        sp[isp].set_ylabel(r"SFR (M$_\odot$/yr)")
        sp[isp].set_yscale('log')
        sp[isp].set_ylim([1.e-5,.2])

    sp[1+nsfr_plots].set_ylabel(r"Wind angle ($^\circ$)")

    sp[0].set_ylabel(r"$\dotM_w$ (M$_\odot$/yr)")
    sp[0].set_yscale('log')

#     sp[3].set_ylabel(r"$\dotM_w/f_\mathrm{edd}$ (M$_\odot$/yr)")
#     sp[3].set_yscale('log')
    
    for this_sp in sp:
#         this_sp.set_xlim([0.,1.2])
        this_sp.set_xlim([0.,0.2])
    
    sp[0].set_ylim([1.e-5,.2])
#     sp[1+nsfr_plots].set_ylim([0.,14.])
    sp[1+nsfr_plots].set_ylim([0.,None])

    sp[nsp-1].set_xlabel("$t$ (Myr)")
    sp[0].legend(loc='best',fontsize='xx-small')
#     sp[0].legend(loc='best')
    
    P.tight_layout()
    P.savefig(outfile)
    P.close()