#import pynbody
import gizmo_tools
import matplotlib.pyplot as plt
import numpy as np

run_id = "3032"

run_names = [
"longrun_medflow_settled_defaultaniso_polar",
"longrun_medflow_vesc_defaultaniso",
"longrun_medflow_vesc_defaultaniso_polar",
"longrun_weakflow_rapid_defaultaniso",
"longrun_weakflow_rapid_defaultaniso_polar",
"longrun_weakflow_settled_defaultaniso",
"longrun_weakflow_settled_defaultaniso_polar",
"longrun_weakflow_vesc_defaultaniso",
"longrun_weakflow_vesc_defaultaniso_polar",
"newflow_settled_thin_up",
"newflow_vesc_thin_45",
"newflow_vesc_thin_side",
"newflow_vesc_thin_up"]

run_names+= ["longrun_medflow_vesc_defaultaniso_polar_highedd"]

run_parameters = gizmo_tools.load_run_parameters("3032")
run_parameters = {x:y for x,y in run_parameters.items() if x in run_names}
ordered_keys = gizmo_tools.run_parameters_names(run_parameters)

run_titles = [x['name'] for x in run_parameters.values()]+["highedd"]

nplots = 4

fig,sp = plt.subplots(nplots,1,figsize=(4,nplots*3),sharex=True,sharey=True)
escfig,escsp = plt.subplots(nplots,1,figsize=(4,nplots*3),sharex=True,sharey=True)

medfig,medsp = plt.subplots()

for run_name,run_title in zip(run_names,run_titles):
    print(f"Loading and plotting {run_name} {run_title}")
    data = np.loadtxt(f"data/t_dyn_frac_{run_id}_{run_name}.dat")
    for iplot in range(nplots):
        sp[iplot].plot(data[:,0]/1.e6,data[:,1+iplot],label=run_title)
        escsp[iplot].plot(data[:,0]/1.e6,data[:,6+iplot],label=run_title)
    medsp.plot(data[:,0]/1.e6,data[:,5],label=run_title)
print("Saving figures")

for iplot,tfactor in enumerate([1,2,4,10]):
    escsp[iplot].set_ylabel(r'fraction $t>{:d}t_{{dyn}}$ or $v_r>v_{{esc}}$'.format(tfactor))
escsp[nplots-1].set_xlabel(r'$t$ (Myr)')
escsp[0].legend(loc='best',fontsize='xx-small')
escfig.savefig("../figures/t_dyn_convergence_vesc.pdf")

for iplot,tfactor in enumerate([1,2,4,10]):
    sp[iplot].set_ylabel(r'fraction $t>{:d}t_{{dyn}}$'.format(tfactor))
sp[nplots-1].set_xlabel(r'$t$ (Myr)')
sp[0].legend(loc='best',fontsize='xx-small')
fig.savefig("../figures/t_dyn_convergence.pdf")

medsp.set_xlabel(r'$t$ (Myr)')
medsp.set_ylabel(r'Median $t/t_{dyn}$')
medsp.set_yscale('log')
medsp.legend(loc='best',fontsize='xx-small')
medfig.savefig("../figures/med_t_dyn.pdf")
