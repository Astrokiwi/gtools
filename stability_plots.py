import gizmo_tools
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


tableDate="180620" 
tableRes="0.1"


boltzmann_cgs = 1.38066e-16

if __name__ == '__main__':

    print("Loading table")
    cloudy_table = gizmo_tools.cloudy_table(tableDate,tableRes,"coolheat_tab_marta/")
    
    print("Setting up input grid")
    L_P = 128
    P_edges = 10**np.linspace(-16,2.,L_P+1)
    P_range = (P_edges[:-1]+P_edges[1:])/2


    L_nH = 128
    nH_edges = 10**np.linspace(0,10,L_nH+1)
    nH_range = (nH_edges[:-1]+nH_edges[1:])/2
    
#     L_T = 32
#     T_edges = 10**np.linspace(0,8,L_T+1)
#     T_range = (T_edges[:-1]+T_edges[1:])/2
    
    n_cells = L_P * L_nH

    table_particles = pd.DataFrame()
    table_particles["AGNIntensity"] = np.full(n_cells,10**1)
    table_particles["AGNDepth"] = np.full(n_cells,5.)
    table_particles["pres"] = np.tile(P_range,L_nH)
    table_particles["nH"] = np.repeat(nH_range,L_P)
    table_particles["temp"] = table_particles["pres"]/table_particles["nH"]/boltzmann_cgs

    print("Calculating dust/cooling/heating properties from table")
    cloudy_table.interp(table_particles)

    table_particles["net_heat"] = (10**table_particles["dHeat"]-10**table_particles["dCool"])/table_particles["nH"]
    
    heating_map = np.reshape(table_particles["net_heat"].data,(L_nH,L_P))
    
    symmetric_max = np.max(np.abs(heating_map))
    
    print("Dumping plot")
    plt.pcolormesh(np.log10(nH_edges),np.log10(P_edges),heating_map.T,cmap='seismic',vmin=-symmetric_max,vmax=symmetric_max,norm=matplotlib.colors.SymLogNorm(1.e-19))
    plt.xlabel('log nH')
    plt.ylabel('log P')
    plt.colorbar()
    plt.savefig("pics/stability.png",dpi=200)
    plt.close('all')