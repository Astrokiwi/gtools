import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P

# for file_base in ["vprof_cold","vprof_allwarm","vprof_warm","vprof_hot","vprof_all"]:

snap_str="100"
for file_base in ["co1","co2","hcn1","hcn2","all"]:

    P.figure()
    prof_data = np.loadtxt("data/line{}{}.dat".format(file_base,snap_str))

    nlines = prof_data.shape[1]-1
    labels = [r"$\phi={}^\circ$".format(x) for x in [0,45,89]]
    for iline in range(nlines):
        P.plot(prof_data[:,0],prof_data[:,iline+1],lw=.5,label=labels[iline])

    P.title("Line: "+file_base+" (angles from edge-on, arbitrary 'flux' units)")
    P.yscale('log')
#     P.ylim([1.e14,None])
    P.xlabel(r"$v_r$ (km/s)")
#     P.ylabel(r"$N_H$ (cm$^{-2}$)")
    P.ylabel(r"$F$ erg/s/MAGICSQUARE/(km/s)")
    P.ylim([1.e46,None])
    P.legend()

    P.savefig("../figures/{}_{}.pdf".format(file_base,snap_str))
    P.close()