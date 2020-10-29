import numpy as np
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as P

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def plot(run_ids,runs,outfile,names=None,sfr_split=None,sfr_plot=True,landscape=False,tsnapshot=None,snap_labels=None,snap_values=None,end_labels=False):
    if names is None:
        names = runs
    if sfr_plot:
        if sfr_split is not None:
            nsfr_plots = len(sfr_split)
        else:
            nsfr_plots = 1
    else:
        nsfr_plots = 0
    if snap_labels is not None and tsnapshot is not None and snap_values is not None:
        nsnap_plots = 1
    else:
        nsnap_plots = 0


    nsp = 2+nsfr_plots+nsnap_plots
    if nsnap_plots>=1:
        nsp+=1 # padding
    
    if type(runs[0])!=list:
        print("just one set")
        runs = [runs]
        names = [names]
    else:
        print("multilist")
    

    if type(run_ids)!=list:
        run_ids = [run_ids for x in runs]
    if type(run_ids[0])!=list:
        run_ids = [[run_ids[i] for x in runs[i]] for i in range(len(runs))]
    else:
        print("individual run_ids per run")
    
#     nsp = 5
#     nsp = 3

#     fig,sp = P.subplots(nsp,1,figsize=(5,8.),sharex=True)
#     fig,sp = P.subplots(nsp,1,figsize=(5,5.),sharex=True)
    ncols = len(run_ids)
    
    scale = 1.
    if landscape:
        fig,allsp = P.subplots(ncols,nsp,figsize=(4.*scale,2.5*ncols*scale),sharex=True,sharey='col') # res version
        allsp = allsp.reshape((ncols,nsp)).T
    else:    
#         fig,allsp = P.subplots(nsp,ncols,figsize=(3.*ncols*scale,4.*scale))#,sharex=True,sharey='row')#) # big paper one
#         height_ratios = []
        if nsnap_plots==0:
            fig,allsp = P.subplots(nsp,ncols,figsize=(3.*ncols*scale,2.*scale*nsp),sharey='row')#,sharex=True,)#) # big paper one
            allsp = allsp.reshape((nsp,ncols))
        else:
            height_ratios = [10]*(nsp-nsnap_plots-1)+[4]+[10]*nsnap_plots
            fig,allsp = P.subplots(nsp,ncols,figsize=(3.*ncols*scale,1.5*scale*nsp),sharey='row',gridspec_kw={'height_ratios':height_ratios})#,sharex=True,)#) # big paper one
            allsp = allsp.reshape((nsp,ncols))
           
    
    skiprate = 1
    
    linewidth = 1.#5 #1.5 #.5

    for icol,(run_ids_set,runs_set,names_set) in enumerate(zip(run_ids,runs,names)):
        sp = allsp[:,icol]
        plot_wind_rate_comparison = False
        if nsnap_plots>0:
            snap_label=snap_labels[icol]
            snap_value=snap_values[icol]
            if snap_label is None or snap_value is None:
                sp[2+nsfr_plots+1].remove()
            else:
                wind_rates_compared = []
                plot_wind_rate_comparison = True
        for irun,(run_id,output_dir,run_label) in enumerate(zip(run_ids_set,runs_set,names_set)):
            print("Plotting: ",output_dir)
        
            if sfr_plot:
                sf_data = np.loadtxt("data/sf"+run_id+output_dir+".dat")
                sf_sp = sp[1]
                if sfr_split is not None:
                    for isplit,sfr_plotlist in enumerate(sfr_split):
                        if run_label in sfr_plotlist:
                            sf_sp  = sp[isplit+1]
                sf_sp.plot(sf_data[:,0][::skiprate]*1.e3,sf_data[:,2][::skiprate],'C'+str(irun),lw=linewidth)
        
            wind_data = np.loadtxt("data/windangle_evolution"+run_id+output_dir+".dat")
            good_slice = (wind_data[:,4]>0.)
            sp[1+nsfr_plots].plot(wind_data[good_slice,1][::skiprate]*1.e3,wind_data[good_slice,2][::skiprate],lw=linewidth)
            skiprate = 1
            p=sp[0].plot(wind_data[good_slice,1][::skiprate]*1.e3,wind_data[good_slice,4][::skiprate],label=run_label,lw=linewidth)
            if end_labels:
                x_l = wind_data[good_slice,1][-1]*1.e3
                y_l = wind_data[good_slice,4][-1]
                sp[0].text(x_l,y_l,run_label,fontsize=6,color=p[0].get_color(),rotation=-10.,rotation_mode='anchor')
            skiprate = 1
    #         sp[3].plot(wind_data[:,1],wind_data[:,4]/edd[irun])
            if plot_wind_rate_comparison:
                rate_index = np.searchsorted(wind_data[good_slice,1]*1.e6,tsnapshot)
                wind_rates_compared.append(wind_data[good_slice,4][rate_index])
        if plot_wind_rate_comparison:
            ncolors = len(snap_value)
            x = np.log10(snap_value)
            y = wind_rates_compared
            sp[2+nsfr_plots+1].scatter(x,y,c=mpl.rcParams["axes.prop_cycle"].by_key()['color'][:ncolors])
            fit = np.polyfit(x, np.log10(y), 1)
            fit_fn = np.poly1d(fit)
            func_label = r"$y=%.2gx^{%.2g}$"%(10.**fit_fn[0],fit_fn[1])
            sp[2+nsfr_plots+1].plot(x,10.**fit_fn(x),ls='--',label=func_label)
#             print(snap_value,wind_rates_compared)
            sp[2+nsfr_plots+1].set_xlabel(snap_label)
            if not end_labels:
                sp[2+nsfr_plots+1].legend(loc='best',fontsize='xx-small')
        if nsnap_plots>=1:
            sp[2+nsfr_plots].set_axis_off() # padding
    
        for isp in range(1,1+nsfr_plots):
            if icol==0:
                sp[isp].set_ylabel(r"SFR (M$_\odot$/yr)")
            sp[isp].set_yscale('log')
            sp[isp].set_ylim([1.e-5,.2])

        if icol==0:
            sp[1+nsfr_plots].set_ylabel(r"Wind angle ($^\circ$)")
            sp[0].set_ylabel(r"$\dotM_w$ (M$_\odot$/yr)")
            if plot_wind_rate_comparison:
                sp[2+nsfr_plots+1].set_ylabel(r"$\dotM_w$ (M$_\odot$/yr)")
#             sp[0].set_ylim([0.,0.01]) # paper
#             sp[0].set_ylim([0.,None]) # paper
#             sp[0].set_ylim([0.,0.001]) # res
#         else:
#             sp[0].set_ylim([0.,0.01]) # paper
#             sp[0].set_ylim([0.,None])
#             sp[0].set_ylim([0.,0.001]) # res
#             if not landscape:
#                 for isp in range(nsp):
#                     sp[isp].yaxis.set_visible(False) # sharey anyway, all good
        sp[0].set_yscale('log')
#         sp[0].set_ylim([1.e-4,0.1])
#         sp[0].set_ylim([1.e-5,1.e-3])
        if plot_wind_rate_comparison:
            sp[2+nsfr_plots+1].set_yscale('log')
#             sp[2+nsfr_plots].set_xscale('log')
            sp[2+nsfr_plots+1].set_ylim([1.e-4,0.1])
        if not landscape:
            for isp in range(0,nsfr_plots+1):
                sp[isp].set_xticklabels("")                    
#                 sp[isp].xaxis.set_visible(False)

    #     sp[3].set_ylabel(r"$\dotM_w/f_\mathrm{edd}$ (M$_\odot$/yr)")
    #     sp[3].set_yscale('log')

        for this_sp in sp[:-1]:
            this_sp.set_xlim([0.,1.e3]) # large scale runs

    
#         for this_sp in sp[:-1]:
    #         this_sp.set_xlim([0.,1.2])
    #         this_sp.set_xlim([0.,0.2])
#             this_sp.set_xlim([0.,0.0099])
#             this_sp.set_xlim([0.,0.0099*1.e3]) # paper
#             this_sp.set_xlim([0.,4.]) # res
    
    #     sp[0].set_ylim([1.e-5,.2])
    #     sp[0].set_ylim([0.,.001])
#         sp[0].set_ylim([0.,None])
#         sp[1].set_ylim([5.,25.])
    #     sp[1+nsfr_plots].set_ylim([0.,14.])

    #     sp[1+nsfr_plots].set_ylim([0.,16.])

        sp[nsp-1-nsnap_plots-1].set_xlabel("$t$ (kyr)")
        sp[nsp-nsnap_plots-1].set_xlabel("$t$ (kyr)")
#         sp[nsp-1].set_xlabel("$t$ (Myr)")
        if not end_labels:
            sp[0].legend(loc='best',fontsize='xx-small')
    #     sp[0].legend(loc='best')
    
#     P.subplots_adjust(hspace=0.,wspace=0.)
    P.tight_layout(h_pad=0.,pad=0.,w_pad=0.)#pad=0.,w_pad=0.1)
    P.savefig(outfile)
    P.close()
    
    
if __name__ == '__main__':
    sfr_split = None

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

#     names = ["02","04"]#,"06","08","1"]
#     runs = ["restest0m"+x for x in names]
#     run_ids     = ["2020"]*3
#     runs = ["restest0m005_tiny","restest0m01_tiny","restest0m02_tiny"]
#     runs = ["restest0m01_small","restest0m02_small","restest0m02","restest0m04"]
#     run_ids = ["2020"]*len(runs)
#     run_ids = ["2020"]*5
#     runs = ["aniso_"+x for x in ["0","1","2","3","4"]]
# #     names = runs
#     outfile = "../figures/timeseries_summary_res.pdf"
#     plot(run_ids,runs,outfile)



#     runs = ['longrun_weakflow_settled_defaultaniso',
#      'longrun_weakflow_vesc_defaultaniso',
#      'longrun_weakflow_rapid_defaultaniso',
#      'longrun_weakflow_settled_defaultaniso_polar',
#      'longrun_weakflow_vesc_defaultaniso_polar',
#      'longrun_weakflow_rapid_defaultaniso_polar',
#      'longrun_medflow_settled_defaultaniso_polar',
#      'longrun_medflow_vesc_defaultaniso',
#      'longrun_medflow_vesc_defaultaniso_polar',
#      'newflow_settled_thin_up',
#      'newflow_vesc_thin_side',
#      'newflow_vesc_thin_45',
#      'newflow_vesc_thin_up']
    runs = [  ["longrun_weakflow_settled_defaultaniso","longrun_weakflow_vesc_defaultaniso","longrun_weakflow_rapid_defaultaniso"],
                ["longrun_weakflow_settled_defaultaniso_polar","longrun_weakflow_vesc_defaultaniso_polar","longrun_weakflow_rapid_defaultaniso_polar"],
                ["longrun_medflow_settled_defaultaniso_polar","longrun_medflow_vesc_defaultaniso","longrun_medflow_vesc_defaultaniso_polar"],
                ["newflow_settled_thin_up","newflow_vesc_thin_side","newflow_vesc_thin_45","newflow_vesc_thin_up"]
                ]


    run_ids = ["3032"]*len(runs)
    outfile = "../figures/timeseries_summary_flows.pdf"
    plot(run_ids,runs,outfile,sfr_plot=False)
    
    