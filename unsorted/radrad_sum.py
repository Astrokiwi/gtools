import numpy as np

for dump in ["000","001"]:
    d = np.loadtxt("data/allp2014q2edd20redo_{}".format(dump),usecols=23)
    radsum = np.sum(d)
    print(dump,radsum)