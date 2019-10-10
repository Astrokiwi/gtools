print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import itertools

import config3032
config3032.setup()

import gizmo_tools

if __name__=='__main__':
    run_id = "3032"
    run_name = "newflow_vesc_thin_45"
    
    line_names_ordered = [r"$F_{IR}$"]+gizmo_tools.lineNames_list
    
    
    print("Run code:",config3032.run_parameters[run_name]["name"])

    lum_codes = ["luminosity"]+["lum_"+line for line in gizmo_tools.line_bases]
    
    fig,sp = plt.subplots()
    ax = sp
    
    for line_name,lum_code in zip(line_names_ordered,lum_codes):
        lumrad_data = np.loadtxt("data/lumrads_{}_{}.dat".format(run_name,lum_code))
        ax.plot(lumrad_data[:,0],lumrad_data[:,1]/np.max(lumrad_data[:,1]),label=line_name)

    ax.set_xlabel(r"$r$ (pc)")    
    ax.set_ylabel(r"$F(<r)/F_{tot}$")    
    ax.legend(prop={'size':8})
    
    fig.savefig("../figures/lumrad_3032.pdf")
    plt.close('all')