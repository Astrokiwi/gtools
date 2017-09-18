import numpy as np
import h5py
# from fof import fof
from fof_step import fof_step
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P
import os

infile = "/export/1/djw/gizmos/1032/goodsurf_fat_fixedtable/snapshot_150.hdf5"

f = h5py.File(infile,"r")

xyz = np.array(f["/PartType0/Coordinates"])
xyz*= 1.e3 # to pc

rad = np.sqrt(xyz[:,0]**2+xyz[:,1]**2)
indisc = (rad<1.5) & (-.25<xyz[:,2]) & (xyz[:,2]<.25)

#indisc = indisc & (xyz[:,1]>0.) & (xyz[:,0]>-0.5) & (xyz[:,0]<0.)

xyz = xyz[indisc,:]

#xyz = xyz[:5000,:]

#xyz[:,2] = 0.

# fof_step.setup(xyz,.005)
# 
# ip = 0
# 
# while(fof_step.step(1)):
#     print(ip)
#     P.clf()
#     x=P.scatter(xyz[:,0],xyz[:,1],c=fof_step.grp,marker='x')
#     P.scatter(xyz[ip,0],xyz[ip,1],c=fof_step.grp[ip],marker='o')
#     P.colorbar(x)
#     P.draw()
#     P.savefig("../pics/fof%04d.png"%ip)
#     ip+=1




# fof_step.step(1)
# x=P.scatter(xyz[:,0],xyz[:,1],c=fof_step.grp,marker='x')
# P.scatter(xyz[ip,0],xyz[ip,1],c=fof_step.grp[ip],marker='o')
# P.colorbar(x)
# P.clim(-1,5)
# 
# def on_keyboard(event):
#     global ip
#     print(event.key)
#     if event.key == ' ':
#         ip+=1
#         fof_step.step(1)
# 
#         P.clf()
#         x=P.scatter(xyz[:,0],xyz[:,1],c=fof_step.grp,marker='x')
#         P.scatter(xyz[ip,0],xyz[ip,1],c=fof_step.grp[ip],marker='o')
#         P.colorbar(x)
#         P.clim(-1,5)
#         P.draw()
# 
# P.gcf().canvas.mpl_connect('key_press_event', on_keyboard)
# 
# P.show()

# fof_step.finish()
# 
# os.system("convert ../pics/fof*.png ../../movies/fof.gif")

#xyz[:,2] = 0.

grp = fof_step.fof(xyz,.002,30)

# grp = fof.calc_fof(xyz,.1)
# 
print(grp)

uniq = np.unique(grp)
print(uniq)

ngrp = -1

# for igrp in uniq:
#     grp_slice = (grp==igrp)
# #     grp_size = np.sum(grp_slice)
# #     if ( grp_size<30 ):
# #         grp[grp_slice] = -1
# #     else:
#     grp[grp_slice] = ngrp
#     ngrp+=1

#print(np.unique(grp))

uniq = np.unique(grp)
print(uniq)
print([[igrp,np.sum(grp==igrp)] for igrp in uniq])



# uniq_grps = np.unique(grp)
# 
# # np.savetxt("all",xyz)
# 
# for igrp in uniq_grps:
#     xyz_out = xyz[(grp==igrp),:]
# #     if ( xyz_out.size//3>2 ):
#     print("output{} size={}".format(igrp,xyz_out.size//3))
#     np.savetxt("group{}".format(igrp),xyz_out)

clump_slice = (grp!=-1)
# clump_slice = (grp!=-2)

fig = P.gcf()
fig.set_size_inches(28,20)

P.clf()
x=P.scatter(xyz[clump_slice,0],xyz[clump_slice,1],c=grp[clump_slice],marker='x',cmap='jet')
cb = P.colorbar(x)
P.xlabel("x (pc)")
P.ylabel("y (pc)")
cb.set_label("clump ID")
#P.draw()
# P.xticks(np.arange(-.5,0.,.01))
# P.yticks(np.arange(.4,1.4,.01))
# P.grid()
P.savefig("../pics/fof_face.png")

# P.clf()
# x=P.scatter(xyz[clump_slice,0],xyz[clump_slice,2],c=grp[clump_slice],marker='x',cmap='tab20c')
# P.colorbar(x)
# #P.draw()
# P.savefig("../pics/fof_side.png")
