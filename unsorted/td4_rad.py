import numpy as np

import gizmo_tools

from sys import path,exit
path.append("../src/")
import tab_interp

# import matplotlib.pyplot as plt
import time


from multiprocessing import Pool
from functools import partial

from functools import wraps


# run_id = "3001"
# run_names = ["a2_e01","a2_e02","a2_e05","a2_e1","a2_e2"]
# edds = [0.01,0.02,0.05,0.1,0.2]
# nruns = len(run_names)
# snap_str = "020"
# max_rad = 5.

bigstart = time.time()
table_time = 0.

run_id = "3032"
run_names = [
# "longrun_medflow_settled_defaultaniso",
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
# edds = [0.01,0.02,0.05,0.1,0.2]
# nruns = len(run_names)

# max_rad = 100.
# min_rad = 1.5

lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
line_codes = ["line_"+line for line in lines]

grain_size = 1. # microns
grain_density = 3 # g/cm**3

lum_factors = np.array([0.5,0.9])
output_linesteps = np.arange(0.,1.01,0.01)

grain_specific_emission_cross_section = 3./(4.*grain_size*1.e-6*grain_density) # cm**2/g, = 250000 !

# grain_specific_emission_cross_section = 1000


grain_specific_emission_cross_section_astronomical = grain_specific_emission_cross_section * 0.000208908219 # to pc**2/Msun
print(grain_specific_emission_cross_section,grain_specific_emission_cross_section_astronomical)

tableDate="060319"
tableRes="0.1"
# chTab = tab_interp.CoolHeatTab( ("coolheat_tab_marta/shrunk_table_labels_"+tableDate+"tau.dat"),
#                                 ("coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_tau.dat"),
#                                 ("coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taunodust.dat"),
#                                 ("coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taunodust.dat"),
#                                 ("coolheat_tab_marta/shrunk_table_labels_"+tableDate+"taudense.dat"),
#                                 ("coolheat_tab_marta/shrunk_table_"+tableDate+"_m"+tableRes+"_hsmooth_taudense.dat")
#                                 )
# interpTabVec = np.vectorize(chTab.interpTab)

cloudy_table = gizmo_tools.cloudy_table(tableDate,tableRes)

# fig = plt.figure()

# fig,ax = plt.subplots(1,1,figsize=(6,6))

def clean_file_not_found(f):
    """ Decorate a function to suppress OSErrors caused by file-not-found in numpy.loadtxt
    """
    
    @wraps(f)
    def g(*args,**kwargs):
        try:
            return f(*args,**kwargs)
        except OSError:
            print(" not found - skipping")
    return g    

@clean_file_not_found
def calc_dump_lumrad(run_name,snap_str):
# for irun,run_name in enumerate(["longrun_weakflow_settled_defaultaniso_polar"]):
    print(run_name,snap_str)
# for run_name in [run_names[0]]:
    header,particles = gizmo_tools.load_gizmo_pandas(run_id,run_name,snap_str,
        ["Masses","Coordinates","Velocities","AGNIntensity","AGNDepth","Density","InternalEnergy","SmoothingLength"])
    gizmo_tools.load_calc_reqs(particles,["nH","temp","rad3d"])

    # this part is slow
    print("tabulating values")
#     start = time.time()
#     tabStructs = interpTabVec(  particles["nH"].astype(np.float64),
#                                 particles["temp"].astype(np.float64),
#                                 particles["AGNIntensity"].astype(np.float64),
#                                 particles["AGNDepth"].astype(np.float64))
# 
#     particles["dustTemp"] = np.array(list(map(lambda y: y.dustT, tabStructs)))
#     particles["dg"] = np.array(list(map(lambda y: y.dg, tabStructs)))

    cloudy_table.interp(particles)
#     end = time.time()
#     table_time+=end-start
#     print("table time:",table_time)
    print(particles.temp.min(),particles.temp.max(),particles.dg.median(),particles.dg.mean(),particles.line_co1.median(),particles.line_co1.mean())


    particles["brightness"] = 5.67e-5 * particles["dustT"]**4. * particles["dg"]/np.nanmax(particles["dg"]) # erg/s/cm^2
    particles["emitArea"] = np.minimum(particles["dg"]*particles["Masses"]*grain_specific_emission_cross_section_astronomical,
                                            np.pi*particles["SmoothingLength"]**2)
    particles["luminosity"] = particles["emitArea"]*particles["brightness"]
    gizmo_tools.load_calc_req(particles,"rad3d")
    
    for line in lines:
#         print(line,particles["line_"+line].min(),particles["line_"+line].max())
        particles["lum_"+line] = particles["line_"+line]*particles["Masses"] # n.b. units irrelevant, we are normalising later

#     gizmo_tools.calculate_vrad(particles)
#     vcut=100.
#     wind_slice = particles["vrad"]>vcut
#     disc_slice = particles["vrad"]<=vcut
#     lum_wind = np.sum(particles["luminosity"][wind_slice])
#     lum_disc = np.sum(particles["luminosity"][disc_slice])
#     print(edds[irun],np.log10(lum_wind),np.log10(lum_disc),lum_wind/lum_disc)


#     particles["luminosity"] = 1. # just mass effectively
    
    
#     label = r"$\eta_\mathrm{{Edd}}={:4.2f}$".format(edds[irun]))
#     label = run_name

    particles.sort_values("rad3d",inplace=True)
    
    summary_outp = []
#     full_outp = [particles["rad3d"].values[::output_nskip]]
#     header_text = "rad"
    for brightness_key in ["luminosity"]+["lum_"+line for line in lines]:
#         print(brightness_key)
#         header_text+=" "+brightness_key
        lum_cum = particles[brightness_key].cumsum()
        lum_tot = lum_cum.iloc[-1]
        
        output_indices = lum_cum.searchsorted(output_linesteps*lum_tot)
        np.savetxt("data/lumrads_{}_{}_{}.dat".format(run_name,brightness_key,snap_str),np.array([particles["rad3d"].values[output_indices],lum_cum.values[output_indices]]).T)
        
        
        r_lum_indices = lum_cum.searchsorted(lum_factors*lum_tot)
        r_lums = particles["rad3d"].iloc[r_lum_indices]
        summary_outp+=[r_lums.values]
#     print("Dumping output")
    np.savetxt("data/summary_lumrads_{}_unextinguished_{}.dat".format(run_name,snap_str),summary_outp)
#     np.savetxt("data/lumrads_{}.dat".format(run_name),np.array(full_outp).T,header=header_text)
#         print(brightness_key,r_lums.values,r_lums.iloc[1]/r_lums.iloc[0])

#     lum_rad = particles["luminosity"].cumsum()
# #     cut_rad_index = particles["rad3d"].searchsorted(max_rad)[0] # pandas <0.24.0
#     cut_rad_index = particles["rad3d"].searchsorted(max_rad) # pandas >=0.24.0
#     cut_rad_inner_index = particles["rad3d"].searchsorted(min_rad) # pandas >=0.24.0
# #     lum_cut = lum_rad.iloc[cut_rad_index]
# #     lum_rad/=lum_cut
#     lum_rad-=lum_rad.iloc[cut_rad_inner_index] # ignore emission from inner region
#     lum_sum=lum_rad.iloc[-1]
#     lum_rad/=lum_sum
#     r0 = particles["rad3d"].min()
#     plt.plot(particles["rad3d"]-r0,lum_rad,label=label)
#     ax.plot(particles["rad3d"],lum_rad,label=label)
    
#     outp = np.array([particles["rad3d"]-r0,lum_rad]).T
#     np.savetxt("../data/edd{:4.2f}.txt".format(edds[irun]),outp,header="% r-ri F/Ftot")
    

# plt.xlim([1.e-2,max_rad])
# plt.xlim([1.,None])
# # plt.xlim([0.,max_rad])
# # ax.set_ylim([0.,1.05])
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# # ax.set_xlabel(r"$r-r_i$ (pc)")    
# ax.set_xlabel(r"$r$ (pc)")    
# ax.set_ylabel(r"$F(<r)/F_{tot}$")    
# ax.legend(prop={'size':8})
# # plt.savefig("../figures/lumrad_test.png")
# # plt.savefig("../figures/lumrad_3001.pdf")
# fig.savefig("../figures/lumrad_3032.pdf")
# plt.close('all')
# 
# bigend = time.time()
# 
# print("table time:",table_time)
# print("total time:",bigend-bigstart)


if __name__=='__main__':
# test
#     snap_str = "100"
#     calc_dump_lumrad(run_names[5],snap_str)

# default for paper
    snap_str = "100"
    for irun,run_name in enumerate(run_names):
        calc_dump_lumrad(run_name,snap_str)


#     v_vs_pos_rays_for_pool = partial(v_vs_pos_rays,run_id,run_name,snap_str)
# #     list(map(v_vs_pos_rays_for_pool,angles))
#     with Pool(processes=3) as pool:
#         pool.map(v_vs_pos_rays_for_pool,angles)

# evolution
#     snap_strs = ["{:03d}".format(x) for x in range(0,310,10)] 
#     for run_name in run_names:
#         one_run_rads = partial(calc_dump_lumrad,run_name)
# #         list(map(one_run_rads,snap_strs))
#         with Pool(processes=64) as pool:
#             pool.map(one_run_rads,snap_strs)
# #         for snap_str in snap_strs:
# #             calc_dump_lumrad(run_name,snap_str)



