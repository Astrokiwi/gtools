from unsorted import time_series_plot


def latex_float(f):
    float_str = "{0:.2g}".format(f)
#     float_str = "{0:.0e}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


# run_names = ["aniso_4","aniso_4_bigrad","aniso_4_biggerrad"]
# run_names = ["superhr","hyperhr"]
# run_ids = ['2020']*len(run_names)
# outfile = "../figures/hr_time_series_summary.png"

# run_names = ["highI_close","highI_close_64n","highI_close_quintic64","highI_close_quintic128"]
# run_ids = ['2022','2022','2024','2024']
# outfile = "../figures/sph_time_series_summary.png"

# run_names = [ ["newtable00001","newtable00005","newtable"],["newtable0005","newtable001","newtable005","newtable01"]]
# run_ids = '2022'
# outfile = "../figures/newtable_res_time_series_summary.pdf"
# res = [[1.e-5,5.e-5,1.e-4],[5.e-4,1.e-3,5.e-3,1.e-2]]
# run_labels = [[r"$"+latex_float(x)+"$ M$_\odot$" for x in y] for y in res]
# time_series_plot.plot(run_ids,run_names,outfile,names=run_labels,sfr_plot=False)

# # run_names = ["fixed_h_restest_0m000001","fixed_h_restest_0m00001","fixed_h_restest_0m0001","fixed_h_restest_0m001","fixed_h_restest_0m01"]
# run_names = ["fixed_h_restest_0m000001","fixed_h_restest_0m00001","fixed_h_restest_0m0001","fixed_h_restest_0m0001_medsoft","fixed_h_restest_0m0001_bigsoft","fixed_h_restest_0m001","fixed_h_restest_0m01"]
# run_ids = ['2022']*len(run_names)
# outfile = "../figures/fixed_h_res_time_series_summary.pdf"
# run_labels = run_names#[[r"$"+latex_float(x)+"$ M$_\odot$" for x in y] for y in res]
# # extract_sfr.extract_sfr(run_ids,run_names)
# # quickangle.dump_full_evolution(run_ids,run_names)
# time_series_plot.plot(run_ids,run_names,outfile,names=run_labels,sfr_plot=False)

# # run_names = ["fixed_h_restest_0m000001","fixed_h_restest_0m00001","fixed_h_restest_0m0001","fixed_h_restest_0m001","fixed_h_restest_0m01"]
# run_names = ["fixed_h_restest_0m0001","holetest_10","correctc","correctc_lowdens"]
# run_names = ["correctc_h05_hmiddens","correctc_h05_lowdens","correctc_h05_middens","correctc2_h05_middens","correctc_h05","correctc_hmiddens","correctc_lowdens","correctc_middens","correctc"]
# run_ids = ['2022']*len(run_names)
# outfile = "../figures/hole_c_tests_time_series_summary.pdf"
# run_labels = run_names#[[r"$"+latex_float(x)+"$ M$_\odot$" for x in y] for y in res]
# # extract_sfr.extract_sfr(run_ids,run_names)
# quickangle.dump_full_evolution(run_ids,run_names)
# time_series_plot.plot(run_ids,run_names,outfile,names=run_labels,sfr_plot=False)

# run_names = ["fixed_h_restest_0m0001","holetest_07","holetest_10","holetest_13","holetest_15"]
# run_ids = ['2022']*len(run_names)
# outfile = "../figures/hole_time_series_summary.pdf"
# run_labels = run_names#[[r"$"+latex_float(x)+"$ M$_\odot$" for x in y] for y in res]
# # extract_sfr.extract_sfr(run_ids,run_names)
# # quickangle.dump_full_evolution(run_ids,run_names)
# time_series_plot.plot(run_ids,run_names,outfile,names=run_labels,sfr_plot=False)

# run_names = ["newtable","newtableN1M4","newtableN2M8","newtableN1M4_np48","newtableN2M8_np48"]
# run_ids = ['2022']*5
# outfile = "../figures/newtable_N_time_series_summary.pdf"

# run_names = ["run_a0_e1","run_a1_e1","run_a2_e01","run_a2_e01_T1000","run_a2_e01_T300","run_a2_e01_T30","run_a2_e05","run_a2_e1","run_a2_e2","run_a3_e1","run_a2_e01_SN100ONE","run_a2_e01_SN1000"]
# run_ids = ['2022']*12
# outfile = "../figures/prodruns_time_series_summary.pdf"

# run_names = ["run_a2_e02"]
# run_ids = ['2022']
# extract_sfr.extract_sfr(run_ids,run_names)
# quickangle.dump_full_evolution(run_ids,run_names)
# time_series_plot.plot(run_ids,run_names,outfile,sfr_plot=False)

# run_names = ["run_a0_e1","run_a1_e1","run_a2_e1","run_a3_e1"]
# run_ids = '2022'
# outfile = "../figures/aniso_prodruns_time_series_summary.pdf"
# time_series_plot.plot(run_ids,run_names,outfile,sfr_plot=False)



# run_names = [   "run_a0_e1","run_a1_e1","run_a2_e1","run_a3_e1",
#                 "run_a2_e01","run_a2_e02","run_a2_e05","run_a2_e1","run_a2_e2",
#                 "run_a2_e01_T1000","run_a2_e01_T300","run_a2_e01_T30","run_a2_e01",
#                 "run_a2_e01","run_a2_e01_SN100","run_a2_e01_SN1000"
#                 ]
# run_names = [   "run_a2_e01_SN100","run_a2_e01_SN1000"
#                 ]
# run_labels = [x[4:] for x in run_names]
# run_names = run_labels # truncated names in 3001
# run_ids = ['3001']*len(run_names)#'2022'

# FOR paper 1
# run_names = [   ["run_a0_e1","run_a1_e1","run_a2_e1","run_a3_e1"],
#                 ["run_a2_e01","run_a2_e02","run_a2_e05","run_a2_e1","run_a2_e2"],
#                 ["run_a2_e01_T1000","run_a2_e01_T300","run_a2_e01_T30","run_a2_e01"],
#                 ["run_a2_e01","run_a2_e01_SN100","run_a2_e01_SN1000"]
#                 ]
# run_labels = [[x[4:] for x in y] for y in run_names]
# run_names = run_labels # truncated names in 3001
# run_ids = '3001'
# # run_ids = '2022'
# outfile = "../figures/all_prodruns_time_series_summary.pdf"
# # outfile = "../figures/all_prodruns_time_series_summary_old.pdf"
# # extract_sfr.extract_sfr(run_ids,run_names)
# # quickangle.dump_full_evolution(run_ids,run_names)
# tsnapshot = 2.e3
# snap_labels = [r"$\log_{10}\eta_a$ $t=2$ kyr",r"$\log_{10}\gamma_\mathrm{Edd}$ $t=2$ kyr",None,None]
# snap_values = [ [1.,10.,100.,1000.],
#                 [0.01,0.02,0.05,0.1,0.2],
#                 None,
#                 None]
# time_series_plot.plot(run_ids,run_names,outfile,names=run_labels,sfr_plot=False,tsnapshot=tsnapshot,snap_labels=snap_labels,snap_values=snap_values)

# FOR TORUS2018
# run_names = [   ["run_a0_e1","run_a1_e1","run_a2_e1","run_a3_e1"],
#                 ["run_a2_e01","run_a2_e02","run_a2_e05","run_a2_e1","run_a2_e2"],
#                 ["run_a2_e01","run_a2_e01_SN100","run_a2_e01_SN1000"]
#                 ]
# run_labels = [[x[4:] for x in y] for y in run_names]
# run_names = run_labels # truncated names in 3001
# run_ids = '3001'
# outfile = "../figures/TORUS2018_time_series_summary.pdf"
# tsnapshot = 2.e3
# snap_labels = [r"$\log_{10}\eta_a$ $t=2$ kyr",r"$\log_{10}\gamma_\mathrm{Edd}$ $t=2$ kyr",None,None]
# snap_values = [ [1.,10.,100.,1000.],
#                 [0.01,0.02,0.05,0.1,0.2],
#                 None,
#                 None]
# time_series_plot.plot(run_ids,run_names,outfile,names=run_labels,sfr_plot=False,tsnapshot=tsnapshot,snap_labels=snap_labels,snap_values=snap_values)

run_names = ["a2_e01_0m00001","a2_e01","a2_e01_0m001","a2_e01_0m01","a2_e01_0m1"]
res = [1.e-5,1.e-4,1.e-3,1.e-2,1.e-1]
outfile = "../figures/res_time_series_summary.pdf"
run_ids = ["3001"]*len(run_names)
run_labels = [r"$"+latex_float(x)+"$ M$_\odot$" for x in res]
# extract_sfr.extract_sfr(run_ids,run_names)
# quickangle.dump_full_evolution(run_ids,run_names)
time_series_plot.plot(run_ids, run_names, outfile, names=run_labels, sfr_plot=False, landscape=True)


# run_names = ["run_a2_e01_T1000","run_a2_e01_T300","run_a2_e01_T30","run_a2_e01"]
# run_ids = ['2022']*4
# outfile = "../figures/floor_prodruns_time_series_summary.pdf"
# time_series_plot.plot(run_ids,run_names,outfile,sfr_plot=False)
# 
# run_names = ["run_a2_e01","run_a2_e05","run_a2_e1","run_a2_e2"]
# run_ids = ['2022']*4
# outfile = "../figures/edd_prodruns_time_series_summary.pdf"
# time_series_plot.plot(run_ids,run_names,outfile,sfr_plot=False)
# 
# run_names = ["run_a2_e01","run_a2_e01_SN100","run_a2_e01_SN1000"]
# run_ids = ['2022']*3
# outfile = "../figures/SN_prodruns_time_series_summary.pdf"
# time_series_plot.plot(run_ids,run_names,outfile,sfr_plot=False)
