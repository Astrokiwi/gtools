print("Importing")

# Import libraries to do our magic

from tools import phaseplots

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P

def latex_float(f):
    float_str = "{0:.2g}".format(f)
#     float_str = "{0:.0e}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


# run_names = ["newtable00001","newtable00005","newtable","newtable0005","newtable001","newtable005","newtable01"]
# res = [1.e-5,5.e-5,1.e-4,5.e-4,1.e-3,5.e-3,1.e-2]
# run_ids = ["2022"]*len(run_names)
# # snap_strs = ["039"]*len(run_names)
# snap_strs = ["015"]*len(run_names)
# doruns = [True]*len(run_names)
# titles = run_names

# run_names = ["fixed_h_restest_0m000001","fixed_h_restest_0m00001","fixed_h_restest_0m0001","fixed_h_restest_0m001","fixed_h_restest_0m01"]
# res = [1.e-6,1.e-5,1.e-4,1.e-3,1.e-2]
# run_ids = ["2022"]*len(run_names)
# snap_strs = ["001"]*len(run_names)
# doruns = [True]*len(run_names)
# titles = run_names

run_names = ["a2_e01_0m00001","a2_e01","a2_e01_0m001","a2_e01_0m01","a2_e01_0m1"]
res = [1.e-5,1.e-4,1.e-3,1.e-2,1.e-1]
run_ids = ["3001"]*len(run_names)
snap_strs = ["020"]*len(run_names)
doruns = [True]*len(run_names)
titles = run_names

nruns = len(run_names)

fig,sp = P.subplots(1,nruns,sharex=True,sharey=True,figsize=(12,3))
# fig,sp = P.subplots(1,nruns,sharex=True,sharey=True,figsize=(36,9))


for irun,(run_id,run_name,snap_str,title,dorun) in enumerate(zip(run_ids,run_names,snap_strs,titles,doruns)):

    if not dorun:
        continue

    loadVals = ["nH_p","TK_p"]

    rcut = 80.


    time,values = phaseplots.loadvalues(run_id, run_name, snap_str, loadVals, rcut)

    print("time=",time)

    this_sp=sp[irun]
    this_sp.set_title(r"$"+latex_float(res[irun])+"$ M$_\odot$")
    mappablePlot = phaseplots.plot_phaseplot(this_sp, values, "nH_p", "TK_p", rcut, cmap='Greys')
#     if irun!=nruns-1:
#         this_sp.set_
    if irun!=0:
        this_sp.yaxis.set_visible(False)

#     cbar = P.colorbar(mappablePlot,label=r"$\log M$ (M$_\odot$)",cax=cbar_ax,orientation='horizontal')

fig.tight_layout()
P.savefig("../figures/res_phaseplot_nH_temp_new{}{}{}.png".format(run_id,run_name,snap_str),dpi=200)
# P.savefig("../figures/res_phaseplot_nH_temp_old{}{}{}.png".format(run_id,run_name,snap_str),dpi=200)
