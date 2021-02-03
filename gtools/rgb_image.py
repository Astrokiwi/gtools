# Import libraries to do our magic
import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as P
import matplotlib.patches as mpatches
# from skimage import exposure

from scipy import ndimage


labels = ["CO(3-2)","CO(6-5)","HCN(4-3)","HCN(8-7)",r"H$_2$ (1-0) S(1)",r"H$_2$ (0-0) S(0)",r"H$_2$ (0-0) S(3)",r"$F_{IR}$"]
line_codes = ["co1m","co2m","hcn1m","hcn2m","h2_1m","h2_2m","h2_3m","view"]

line_label_dict = dict(zip(line_codes,labels))

def produce_and_save_rgb_table_image(files,ax=None,outfile=None,rad=None,labels=None,intensities=[1.,1.,1.],max_contrast=False,cmap='plasma',vrange=None,printrange=False,gauss_rad=None): #
    rgb_mode = True
    if len(files)==1:# and len(intensities)==1:
        rgb_mode = False
    else:
        if len(files)!=3:
            raise Exception("Need exactly three files for RGB image")
#         if len(intensities)!=3:
#             raise Exception("Need exactly three intensities if specified")
#     intensities = np.array(intensities)
    if (ax is None and outfile is None) or (ax is not None and outfile is not None):
        raise Exception("need to specify *either* axis or outfile")
    vmin = None
    vmax = None
    if vrange is not None:
        vmin = vrange[0]
        vmax = vrange[1]
    maps = [np.loadtxt(file) for file in files]

    gauss_pix = None
    if rad is None:
        extent = None
        if gauss_rad is not None:
            gauss_pix = gauss_rad
    else:
        extent = [-rad,rad,-rad,rad]
        if gauss_rad is not None:
            gauss_pix = int(gauss_rad/(2*rad)*maps[0].shape[0]) # assume square

    for imap,map in enumerate(maps):
#         finite_slice = np.isfinite(map)
#         map-=np.nanmin(map[finite_slice])
#         print(np.nanmin(map[finite_slice]),np.nanmax(map[finite_slice]))
#         map/=np.nanmax(map[finite_slice])
#         map=np.arcsinh(map)
#         if imap==2:
#         map = exposure.equalize_hist(map)
        map*=intensities[imap]
        print(np.nanmin(map),np.nanmax(map))
        maps[imap]=map
        if gauss_pix is not None:
            maps[imap] = ndimage.gaussian_filter(maps[imap],gauss_pix)

        
    rgb_map = np.array([map.T for map in maps]).T # to size x size x 3
#     rgb_map = exposure.equalize_hist(rgb_map)
#     if not rgb_mode:
#         rgb_map = np.tile(rgb_map.T,(3,1,1)).T
#     if rgb_mode:
#         for imap in range(3):
#             rgb_map[:,:,imap]*=intensities[imap]

    if ax is None:
        P.clf()
        fig, ax = P.subplots()
    if rgb_mode:
        mappable=ax.imshow(rgb_map,extent=extent,origin='lower')
    else:
        mappable=ax.imshow(maps[0],extent=extent,origin='lower',cmap=cmap,vmin=vmin,vmax=vmax)
    if extent is not None:
        ax.set_xlabel("pc")
        ax.set_ylabel("pc")
    if labels is not None:
        if len(labels)==1:
            ax.set_title(label=labels[0])
        else:
            red_patch = mpatches.Patch(color=(1.,0.,0.), label=labels[0])
            blue_patch = mpatches.Patch(color=(0.,1.,0.), label=labels[1])
            green_patch = mpatches.Patch(color=(0.,0.,1.), label=labels[2])
            ax.legend(handles=[red_patch,blue_patch,green_patch],fontsize='xx-small')

    if outfile is not None:
        P.savefig(outfile,dpi=300)
#     P.savefig(outfile,dpi=100)
    return mappable

def plot_all():
    intensities = [1.,1.,1.]
    rad = 1.
    basefile = "../data/smoothsum_co1mco2mhcn1mhcn2mh2_1mh2_2mh2_3mgiz_3001_a2_e01_000_%03d.dat"
#     for co_id in [1]:
#         for hcn_id in [3]:
#             for h2_id in [5]:
    for co_id in range(2):
        for hcn_id in range(2,4):
            for h2_id in range(4,7):
                files = [   basefile%hcn_id,
                            basefile%co_id,
                            basefile%h2_id]
                outfile = "../../figures/co_%01d_hcn_%01d_h2_%01d_rgb_boosted.png"%(co_id+1,hcn_id-1,h2_id-3)
                plot_labels = [labels[i] for i in [hcn_id,co_id,h2_id]]
                produce_and_save_rgb_table_image(files,outfile=outfile,rad=rad,intensities=intensities,labels=plot_labels,max_contrast=True)

def plot_nice():
    basefile = "../data/smoothsum_co1mco2mhcn1mhcn2mh2_1mh2_2mh2_3mgiz_3001_a2_e01_000_%03d.dat"
    intensities = [1.,1.,1.]
    rad = 1.
    fig,ax=P.subplots(1,4,sharex=True,sharey=True,figsize=(12,4))
    
    toplot = [3,1,4]
    for iax,iscalar_plot in enumerate(toplot):
        label = [labels[iscalar_plot]]
        infile = [basefile%iscalar_plot]
        produce_and_save_rgb_table_image(infile,ax=ax[iax],rad=rad,intensities=intensities,labels=label,max_contrast=True)
        if iax!=0:
            ax[iax].set_ylabel("")
    all_label = [labels[x] for x in toplot]
    infiles = [basefile%x for x in toplot]
    produce_and_save_rgb_table_image(infiles,ax=ax[3],rad=rad,intensities=intensities,labels=all_label,max_contrast=True)
    ax[3].set_ylabel("")
#     P.subplots_adjust(wspace=0.1)
    P.savefig("../../figures/lines_together.png",dpi=200)

if __name__ == '__main__':
    # Data produced with:
    # python anim.py 3001 a2_e01 --rad 1. --savemap --view face --phi 90. --snap0 100 --maxsnapf 100 --plot co1m,co2m,hcn1m,hcn2m,h2_1m,h2_2m,h2_3m
    
#     plot_all()

    plot_nice()















