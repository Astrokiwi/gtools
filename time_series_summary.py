import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P

if __name__ == '__main__':
#     run_id = "2014"
#     run_ids     = ["2014","2014","2015","2014","2015","2014"]
#     runs = ["q2redo","q2edd05redo","q2edd10_aniso1fixed","q2edd10redo","q2edd10_aniso3fixed","q2edd20redo"]
# #     edd = [.01,.05,.1,.1,.1,.2] # incorrect anyway
#     names = ["A","B","C1","C2","C3","D"]

#     run_ids     = ["1055","1055"]
#     runs = ["q2_SN","q2_SN_slow"]
#     names = ["ASS","AS"]
    run_ids     = ["1058","1058"]
    runs = ["restest0m02","restest0m04"]
    names = ["0m02","0m04"]
    
    
    nsp = 3

    fig,sp = P.subplots(nsp,1,figsize=(6,3.*nsp),sharex=True)

    for irun,(run_id,output_dir) in enumerate(zip(run_ids,runs)):
        print("Plotting: ",output_dir)
        
        sf_data = np.loadtxt("data/sf"+run_id+output_dir+".dat")
        sp[1].plot(sf_data[:,0],sf_data[:,2],lw=.7)
        
        wind_data = np.loadtxt("data/windangle_evolution"+run_id+output_dir+".dat")
        good_slice = (wind_data[:,4]>0.)
        sp[2].plot(wind_data[good_slice,1],wind_data[good_slice,2],lw=.7)
        sp[0].plot(wind_data[good_slice,1],wind_data[good_slice,4],label=names[irun],lw=.7)
#         sp[3].plot(wind_data[:,1],wind_data[:,4]/edd[irun])
    
    sp[1].set_ylabel(r"SFR (M$_\odot$/yr)")
    sp[1].set_yscale('log')

    sp[2].set_ylabel(r"Wind angle ($^\circ$)")

    sp[0].set_ylabel(r"$\dotM_w$ (M$_\odot$/yr)")
    sp[0].set_yscale('log')

#     sp[3].set_ylabel(r"$\dotM_w/f_\mathrm{edd}$ (M$_\odot$/yr)")
#     sp[3].set_yscale('log')
    
    for this_sp in sp:
        this_sp.set_xlim([0.,1.2])
    
    sp[0].set_ylim([1.e-5,.2])
    sp[1].set_ylim([1.e-5,.2])
    sp[2].set_ylim([0.,14.])

    sp[nsp-1].set_xlabel("$t$ (Myr)")
#     sp[0].legend(loc='best',fontsize='xx-small')
    sp[0].legend(loc='best')
    
    P.tight_layout()
#     P.savefig("../figures/timeseries_summary_q2.pdf")
#     P.savefig("../figures/timeseries_summary_SN.pdf")
    P.savefig("../figures/timeseries_summary_res_fast.pdf")
    P.close()