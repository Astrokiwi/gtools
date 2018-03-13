import numpy as np
import calendar
import datetime
import os

now = datetime.datetime.now()

month = now.month
year = now.year
day = now.day

# lastmonth = month-1
# if ( lastmonth == 0 ):
#     lastmonth = 12
#     lastmonthyear-=1
# else:
#     lastmonthyear = year
# 
# days_lastmonth = calendar.monthrange(lastmonthyear,lastmonth)[1]
#sar -f /var/log/sa/sa29

# just do this month for now

files = " "

for iday in range(1,day+1):
    os.system("sar -f /var/log/sa/sa{0:02d} > data/sar_out{0:02d}".format(iday))
    files = files + " data/sar_out{0:02d}".format(iday)
    
sumfile = "data/sar_cat_{:02d}{:02d}.dat".format(month,year)

os.system("cat {} > {}".format(files,sumfile))

allsar_d = np.genfromtxt(sumfile, usecols = (2,3,4,5,6,7), invalid_raise=False)

notnan = ~np.any(np.isnan(allsar_d),axis=1)
allsar_d = allsar_d[notnan,:]

nvalues = allsar_d.shape[0]

timeago = np.flipud(np.arange(0,nvalues*10.,10.))

output = np.vstack([timeago,allsar_d.T])

outputfile = "data/sar_sum_{:02d}{:02d}.dat".format(month,year)


np.savetxt(outputfile,output.T,fmt="%5.2f",header="MinutesAgo %user %nice %system %iowait %steal %idle")