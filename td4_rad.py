import numpy as np

import gizmo_tools

from sys import path
path.append("src/")
import tab_interp

import matplotlib.pyplot as plt

run_id = "3001"
run_names = ["a2_e01","a2_e02","a2_e05","a2_e1","a2_e2"]
# run_names = ["a2_e01","a2_e02","a2_e05"]
edds = [0.01,0.02,0.05,0.1,0.2]
nruns = len(run_names)

# max_rad = 40.
max_rad = 5.

grain_size = 1. # microns
grain_density = 3 # g/cm**3



grain_specific_emission_cross_section = 3./(4.*grain_size*1.e-6*grain_density) # cm**2/g, = 250000 !

# grain_specific_emission_cross_section = 1000


grain_specific_emission_cross_section_astronomical = grain_specific_emission_cross_section * 0.000208908219 # to pc**2/Msun
print(grain_specific_emission_cross_section,grain_specific_emission_cross_section_astronomical)

tableDate="281118"
tableRes="0.0001"
chTab = tab_interp.CoolHeatTab( ("coolheat_tab_marta/shrunk_table_labels_"+tableDate+"tau.dat"),
                                ("coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_tau.dat"),
                                ("coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taunodust.dat"),
                                ("coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taunodust.dat"),
                                ("coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taudense.dat"),
                                ("coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taudense.dat")
                                )
interpTabVec = np.vectorize(chTab.interpTab)

# fig = plt.figure()

for irun,run_name in enumerate(run_names):
#     print(irun,run_name)
# for run_name in [run_names[0]]:
#     snap_str = "100"
    snap_str = "020"
    header,particles = gizmo_tools.load_gizmo_pandas(run_id,run_name,snap_str,
        ["Masses","Coordinates","Velocities","AGNIntensity","AGNDepth","Density","InternalEnergy","SmoothingLength"])
    gizmo_tools.load_calc_reqs(particles,["nH","temp","rad3d"])

    tabStructs = interpTabVec(  particles["nH"].astype(np.float64),
                                particles["temp"].astype(np.float64),
                                particles["AGNIntensity"].astype(np.float64),
                                particles["AGNDepth"].astype(np.float64))

    particles["dustTemp"] = np.array(list(map(lambda y: y.dustT, tabStructs)))
    particles["dg"] = np.array(list(map(lambda y: y.dg, tabStructs)))
    particles["brightness"] = 5.67e-5 * particles["dustTemp"]**4. * particles["dg"]/np.nanmax(particles["dg"]) # erg/s/cm^2
    particles["emitArea"] = np.minimum(particles["dg"]*particles["Masses"]*grain_specific_emission_cross_section_astronomical,
                                            np.pi*particles["SmoothingLength"]**2)
    particles["luminosity"] = particles["emitArea"]*particles["brightness"]
    gizmo_tools.calculate_vrad(particles)
#     vcut=100.
#     wind_slice = particles["vrad"]>vcut
#     disc_slice = particles["vrad"]<=vcut
#     lum_wind = np.sum(particles["luminosity"][wind_slice])
#     lum_disc = np.sum(particles["luminosity"][disc_slice])
#     print(edds[irun],np.log10(lum_wind),np.log10(lum_disc),lum_wind/lum_disc)


#     particles["luminosity"] = 1. # just mass effectively
    
    
    particles.sort_values("rad3d",inplace=True)

    lum_rad = particles["luminosity"].cumsum()
    cut_rad_index = particles["rad3d"].searchsorted(max_rad)[0]
#     lum_cut = lum_rad.iloc[cut_rad_index]
#     lum_rad/=lum_cut
    lum_sum = lum_rad.iloc[-1]
    lum_rad/=lum_sum
    r0 = particles["rad3d"].min()
    plt.plot(particles["rad3d"]-r0,lum_rad,label=r"$\eta_\mathrm{{Edd}}={:4.2f}$".format(edds[irun]))
    
    outp = np.array([particles["rad3d"]-r0,lum_rad]).T
    np.savetxt("../data/edd{:4.2f}.txt".format(edds[irun]),outp,header="% r-ri F/Ftot")
    

# plt.xlim([1.e-2,max_rad])
plt.xlim([0.,max_rad])
plt.ylim([0.,1.05])
# plt.xscale('log')
# plt.yscale('log')
plt.xlabel(r"$r-r_i$ (pc)")    
plt.ylabel(r"$F(<r)/F_{tot}$")    
plt.legend()
# plt.savefig("../figures/lumrad_test.png")
plt.savefig("../figures/lumrad_3001.pdf")
plt.close('all')
