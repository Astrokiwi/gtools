import numpy as np
import matplotlib.pyplot as plt

run_name = "longrun_medflow_vesc_defaultaniso_polar"
run_id = "3032"
snap_str = "200"
val = "nH_p"

# lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
# line_labels = [r"CO(3-2), $866.727$ $\mu$m","CO(6-5), $433.438$ $\mu$m","HCN(4-3), $845.428$ $\mu$m","HCN(8-7), $422.796$ $\mu$m",r"H$_2$ (1-0) S(1), $2.121$ $\mu$m",r"H$_2$ (0-0) S(0), $28.18$ $\mu$m",r"H$_2$ (0-0) S(3), $9.66$ $\mu$m"]
lines = ["co1","co2","hcn1","hcn2","h2_1","h2_3","h2_2"]
line_labels = [r"CO(3-2), $866.727$ $\mu$m","CO(6-5), $433.438$ $\mu$m","HCN(4-3), $845.428$ $\mu$m","HCN(8-7), $422.796$ $\mu$m",r"H$_2$ (1-0) S(1), $2.121$ $\mu$m",r"H$_2$ (0-0) S(3), $9.66$ $\mu$m",r"H$_2$ (0-0) S(0), $28.18$ $\mu$m"]

nlines = len(lines)
ncuts = 3
cut_titles = [r"$T_d<15$ K",r"$15\leq T_d<100$ K",r"$T_d\geq100$ K"]

scale = .7
fig,sp = plt.subplots(ncuts,1,sharex=True,sharey=True,figsize=(6*scale,ncuts*3*scale),constrained_layout=True)

for icut,cut_title in enumerate(cut_titles):
    for iline,(line,line_label) in enumerate(zip(lines,line_labels)):
        infile_name = f"data/median_emmissivity_{run_id}_{run_name}_{snap_str}_{val}_{line}_{icut}.dat"
        x,y = np.loadtxt(infile_name,unpack=True)
        sp[icut].plot(x,y-np.nanmax(y),label=line_label)
        sp[icut].set_title(cut_title)
        sp[icut].set_ylabel(r"$\log_{10}(\epsilon/\epsilon_{max})$ (erg g$^{-1}$ s$^{-1}$)")
sp[2].legend(loc='best',fontsize='x-small')
sp[2].set_xlabel(r"$\log_{10}n_H$ (cm$^{-3}$)")

fig.savefig("../figures/line_behaviour.pdf")