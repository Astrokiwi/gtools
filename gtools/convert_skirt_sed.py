#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

from astropy import constants
from astropy import units

def ryd2microns(x):
    return (1./(constants.Ryd*x)).to(units.um).value

if __name__ == '__main__' :
    infile = sys.argv[1]
    outfile = sys.argv[2]

    nu,nuInu = np.loadtxt(infile,usecols=[0,1],unpack=True)

    df = pd.DataFrame()
    df["wlength"] = ryd2microns(nu)
    df["I_wlength"] = nuInu/df["wlength"]

    df.sort_values(inplace=True,by="wlength")

    df.to_csv(outfile,sep=" ",header=None,index=False)