import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P

if __name__ == '__main__':
    run_id = "2014"
    runs = ["q2redo","q2edd05redo","q2edd10_aniso1redo","q2edd10redo","q2edd10_aniso3redo","q2edd20redo"]
    edd = [.01,.05,.1,.1,.1,.2]
    names = ["A","B","C1","C2","C3","D"]
    
    
    nsp = 4

    fig,sp = P.subplots(nsp,1,figsize=(4,2.*nsp),sharex=True)

    for irun,output_dir in enumerate(runs):
        print("Plotting: ",output_dir)
        sf_data = np.loadtxt("data/sf"+run_id+output_dir+".dat")
        sp[0].plot(sf_data[:,0],sf_data[:,2])
        
        wind_data = np.loadtxt("data/windangle_evolution"+run_id+output_dir+".dat")
        sp[1].plot(wind_data[:,1],wind_data[:,2])
        sp[2].plot(wind_data[:,1],wind_data[:,4],label=names[irun])
        sp[3].plot(wind_data[:,1],wind_data[:,4]/edd[irun])
    
    sp[0].set_ylabel(r"SFR (M$_\odot$/yr)")
    sp[0].set_yscale('log')

    sp[1].set_ylabel(r"Wind angle ($^\circ$)")

    sp[2].set_ylabel(r"$\dotM_w$ (M$_\odot$/yr)")
    sp[2].set_yscale('log')

    sp[3].set_ylabel(r"$\dotM_w/f_\mathrm{edd}$ (M$_\odot$/yr)")
    sp[3].set_yscale('log')

    sp[nsp-1].set_xlabel("$t$ (Myr)")
    sp[2].legend(loc='best',fontsize='xx-small')
    
    P.tight_layout()
    P.savefig("../figures/timeseries_summary_q2.pdf")
    P.close()