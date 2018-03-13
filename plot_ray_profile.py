import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P

for file_base in ["vprof_cold","vprof_allwarm","vprof_warm","vprof_hot","vprof_all"]:

    P.figure()
    prof_data = np.loadtxt("data/{}.dat".format(file_base))

    nlines = prof_data.shape[1]-1
    for iline in range(nlines):
        P.plot(prof_data[:,0],prof_data[:,iline+1],lw=.5)

    P.yscale('log')
    P.ylim([1.e14,None])
    P.xlabel(r"$v_r$ (km/s)")
    P.ylabel(r"$N_H$ (cm$^{-2}$)")

    P.savefig("../figures/{}.pdf".format(file_base))
    P.close()