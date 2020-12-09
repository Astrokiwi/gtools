import tblib.pickling_support
tblib.pickling_support.install()

# Import libraries to do our magic
import numpy as np
import os
this_dir, this_filename = os.path.split(__file__)

# from joblib import Parallel, delayed
from multiprocessing import Pool

from . import frame

from . import gizmo_tools

import argparse

import tqdm

# def sort_nicely( l ):
#     """ Sort the given list in the way that humans expect.
#     """
#     convert = lambda text: int(text) if text.isdigit() else text
#     alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
#     l.sort( key=alphanum_key )

def dump_rad0(infile):
    import h5py
    with h5py.File(infile,"r") as f:
        xyz_p = np.array(f["/PartType0/Coordinates"]) # kpc
        id_p = np.array(f["/PartType0/ParticleIDs"]).astype(int)
    id_p-=1
    rad_p = np.sqrt(xyz_p[:,0]**2+xyz_p[:,1]**2+xyz_p[:,2]**2)*1.e3
    rad_out = np.zeros_like(rad_p)
    rad_out[id_p] = rad_p
    np.save("rad0.npy",rad_out)
#     sys.exit()
    
def string_to_list_or_float(s):
    if s is None:
        return None
    if "," in s:
        return [float(x) for x in s.split(",")]
    return float(s)

# joblib doesn't dump exceptions well, just print them
def plotter_parallel_exception_wrapper(*args,**kwargs):
    try:
        return frame.makesph_trhoz_frame(*args, **kwargs)
    except Exception as e:
        print(e) # should give some output at least
        raise(e) # might not work

def parse_input():
    default_values = dict()
    default_values["nprocs"]=8
    default_values["maxsnapf"]=-1
    default_values["snap0"]=-1
    default_values["rad"]="15"
    #default_values["rads"]=""
    default_values["L"]=400
    default_values["plot"]="dens,temp"
    default_values["views"]="face,side"
    default_values["cmap"]="viridis"
    default_values["slice"]=False
    default_values["phi"]="0."
    default_values["theta"]=0.
    default_values["noaxes"]=False
    default_values["centredens"]=False
    default_values["centrecom"]=False
    default_values["suffix"]=""
    default_values["dotmode"]=""
    default_values["absurd"]=False
    default_values["pixsize"]=1
    default_values["noring"]=False
    default_values["data_ranges"]=""
    default_values["savemap"]=False
    default_values["gaussian"]=None
    default_values["gizmoDir"]=None
    default_values["snapstep"]=1
    default_values["verbose"]=False
    default_values["dataDir"] = None
#    default_values["opac_mu"]=None
    
    parsevals = ["data_ranges","nprocs","maxsnapf","run_id","output_dir","plot","cmap","rad","L","slice","views","phi","theta","noaxes","centredens","centrecom","suffix","dotmode","absurd","pixsize","noring","snap0","savemap","gaussian","gizmoDir","snapstep","verbose","dataDir"]
    #,"opac_mu"]

    parser = argparse.ArgumentParser()
    parser.add_argument('run_id',help="name of superdirectory for runs")
    parser.add_argument('output_dir',help="name of subdirectory for run")
    parser.add_argument('--nprocs',type=int,help="processors to run on (default {})".format(default_values["nprocs"]))
    parser.add_argument('--maxsnapf',type=int,help="snapshot to end on (default=-1=do all snapshots)")
    parser.add_argument('--snap0',type=int,help="snapshot to start on (default=-1=do all snapshots)")
    parser.add_argument('--snapstep',type=str,help="snapshop step rate")
    parser.add_argument('--rad',type=str,help="radius of all plots in parsecs - same value for all plots, or separated by commas")
    #parser.add_argument('--rads',type=float,help="radius of each plot in parsecs, separated by commas")
    parser.add_argument('--L',type=int,help="size of plot area in pixels")
    parser.add_argument('--plot',type=str,help="values to plot, separated by commas")
    parser.add_argument('--views',type=str,help="face and/or side view, separated by commas")
    parser.add_argument('--cmap',type=str,help="colourmap palette")
    parser.add_argument('--phi',type=str,help="phi (latitudinal) rotation angle in degrees")
    parser.add_argument('--theta',type=float,help="theta (azimuthal) rotation angle in degrees")
    parser.add_argument('--slice',help="option to be slice plot",action='store_true')
    parser.add_argument('--noaxes',help="don't show axes - fill frame with plot",action='store_true')
    parser.add_argument('--centredens',help="centre on densest particle",action='store_true')
    parser.add_argument('--centrecom',help="centre on centre of mass",action='store_true')
    parser.add_argument('--suffix',help="suffix on anim filename")
    parser.add_argument('--dotmode',help="max, min, or nothing - dot plot mode")
    parser.add_argument('--pixsize',type=int,help="pixelation level")
    parser.add_argument('--data_ranges',type=str,help="explicit data bounds for plot - format: --data_ranges=min1,max1+min2,max2+min3,max3 etc; the equals sign may be necessary!")
    parser.add_argument('--absurd',help="I'll try spinning, that's a good trick",action='store_true')
    parser.add_argument('--noring',help="Edge on ring, not wrapped",action='store_true')
    parser.add_argument('--savemap',help="Save maps as data file",action='store_true')
    parser.add_argument('--gaussian',type=str,help="Gaussian filter size - same value for all plots, or separated by commas")
    parser.add_argument('--gizmoDir',type=str,help="Custom directory for gizmo runs")
    parser.add_argument('--verbose',help="Print more text",action='store_true')
    parser.add_argument('--dataDir', type=str, help="Custom directory for output data")
#    parser.add_argument('--opac_mu',type=float,help="Wavelength in microns of dust opacity to use")
    args = parser.parse_args()
    
    anim_prams = dict()
    
    # copy from parser into run_prams dict
    for parseval in parsevals:
        if ( vars(args)[parseval] ):
            anim_prams[parseval] = vars(args)[parseval]
            if args.verbose:
                print("setting {} to {}, current value:{}".format(parseval,vars(args)[parseval],anim_prams[parseval]))
        else:
            if ( parseval in default_values ):
                anim_prams[parseval] = default_values[parseval]
                if args.verbose:
                    print("setting {} to default {}, current value:{}".format(parseval,default_values[parseval],anim_prams[parseval]))
            else:
                raise Exception("No default value for {} - it must be specified!".format(parseval))
    return anim_prams


# TODO: copy other keys into frame_prams kwargs, to simplify function call
# and/or use functools properly
def frame_i(i):
    global infiles,outfiles,rads,thetas,phis,frame_prams
#         return sph_frame.makesph_trhoz_frame(infiles[i],outfiles[i],cmap=cmap,flat=flatPlot,ring=ringPlot,plot=toplot,L=L,scale=rads[i],views=toview
    #,rot=[thetas[i],phis[i]],visibleAxes=visibleAxes,centredens=centredens,centrecom=centrecom,dotmode=dotmode,pixsize=pixsize
    #,data_ranges=data_ranges,return_maps=savemap,gaussian=gaussian)
    return frame.makesph_trhoz_frame(infiles[i], outfiles[i], **frame_prams, scale=rads[i], rot=[thetas[i], phis[i]])
    #,opac_mu=opac_mu)


def animate(anim_prams):
    global infiles,outfiles,rads,thetas,phis,frame_prams # I hate this, but it makes pickle happy

    frame_prams = dict()

    frame_prams["flat"] = not anim_prams["slice"]
    frame_prams["ring"] = frame_prams["flat"]
    if anim_prams["noring"]:
        frame_prams["ring"] = False

    if anim_prams["dotmode"]:
        smooth_str = anim_prams["dotmode"]
    elif frame_prams["flat"]:
        smooth_str = "smooth"
    else:
        smooth_str = "slice"
    
    frame_prams["gaussian"] = string_to_list_or_float(anim_prams["gaussian"])
    anim_prams["phi"] = string_to_list_or_float(anim_prams["phi"])

    frame_prams["return_maps"] = anim_prams["savemap"]
    
    if len(anim_prams["data_ranges"])<=0:
        frame_prams["data_ranges"] = None
    else:
        frame_prams["data_ranges"]=[[float(y) for y in x.split(",")] for x in anim_prams["data_ranges"].split("+")]
    
    frame_prams["plot"] = anim_prams["plot"].split(",")
    outp_plot = "".join(frame_prams["plot"])

    frame_prams["views"] = anim_prams["views"].split(",")
#     outp_views = "".join(frame_prams["views"])
    
    if isinstance(anim_prams["phi"],list):
        anim_prams["phi"]=[x*2.*np.pi/360. for x in anim_prams["phi"]]
    else:
        anim_prams["phi"]*=2.*np.pi/360.
    anim_prams["theta"]*=2.*np.pi/360.
    
    frame_prams["visibleAxes"]=not anim_prams["noaxes"]
    
    prams_to_copy = ["cmap","L","centredens","centrecom","dotmode","pixsize"]
    for key in prams_to_copy:
        frame_prams[key] = anim_prams[key]
    
    gizmoDir = anim_prams["gizmoDir"]
    run_id = anim_prams["run_id"]
    output_dir = anim_prams["output_dir"]
    snap0 = anim_prams["snap0"]
    maxsnapf = anim_prams["maxsnapf"]
    snapstep = anim_prams["snapstep"]
    absurd = anim_prams["absurd"]
    rad = string_to_list_or_float(anim_prams["rad"])



    if gizmoDir is None:
        gizmoDir = gizmo_tools.getGizmoDir(run_id)
    movieDir = gizmo_tools.getMovieDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir

    snapi = 0
    if snap0>0:
        snapi=snap0
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir,False,gizmoDir=gizmoDir)

    
# max run
#     if ( maxsnapf>-1 and snapf>maxsnapf ):
    if maxsnapf>-1:
        print("Forcing snapf from {} to {}".format(snapf,maxsnapf))
        snapf = maxsnapf

    if anim_prams["dataDir"] is None:
        data_dir = gizmo_tools.getDataDir()
    else:
        data_dir = anim_prams["dataDir"]
    pic_dir = gizmo_tools.getPicDir()
    os.system("rm "+pic_dir+"/sphplot"+run_id+output_dir+"???.png")

#     print("nfiles:",snapf-snapi+1)
    
    isnaps = range(snapi,snapf+1,snapstep)
    
    infiles = [fullDir+"/snapshot_"+("%03d" % snapx)+".hdf5" for snapx in isnaps]
    outfiles = [pic_dir+"/sphplot"+run_id+output_dir+"%03d.png"%snapx for snapx in isnaps]
    
    nfiles = len(infiles)
    if absurd:
        period2 = 11
        phis = np.linspace(0.,np.pi/period2*nfiles,nfiles)
        thetas = np.linspace(0.,np.pi/period2*nfiles,nfiles)
        period_alt = 10
        ungles = np.linspace(0.,np.pi/period_alt*nfiles,nfiles)
        rads = np.sin(ungles)*rad/2.+rad
    else:
        thetas = np.full(nfiles,anim_prams["theta"])
        if isinstance(anim_prams["phi"],list):
            phis=[anim_prams["phi"]]*nfiles
        else:
            phis = np.full(nfiles,anim_prams["phi"])
        if isinstance(rad,list):
            rads = [rad]*nfiles
        else: # if float
            rads = np.full(nfiles,rad)
    
    if "rad0" in frame_prams["plot"]:
        dump_rad0(infiles[0])

    with Pool(processes=anim_prams["nprocs"]) as pool:
        maps=[]
        for _ in tqdm.tqdm(pool.imap_unordered(frame_i,range(snapf+1-snapi)),total=snapf+1-snapi):
            maps.append(_)

    for result in maps:
        if isinstance(result, frame.ExceptionWrapper):
            result.re_raise()
    if anim_prams["savemap"]:
        plot_strs = list(np.repeat(frame_prams["plot"],len(frame_prams["views"])))

    
        for itime,maps_timeslice in enumerate(maps):
            idump = isnaps[itime]
            #print(len(maps_timeslice))
            for iplot,map in enumerate(maps_timeslice):
                #print(iplot,idump,frame_prams["plot"],plot_strs)
                plot_str = plot_strs[iplot]
                if frame_prams["plot"].count(plot_str)>0:
                    plot_str+="_%03d"%plot_strs[:iplot].count(plot_strs[iplot])
                #print("type=",type(map))
                if type(map) is list:
                    for imap,submap in enumerate(map):
                        outfile = data_dir+"/"+smooth_str+"sum_giz_"+run_id+"_"+output_dir+anim_prams["suffix"]+"_%03d"%idump+"_"+plot_str+"_%03d.dat"%imap
                        np.savetxt(outfile,submap)
                else:
                    outfile = data_dir+"/"+smooth_str+"sum_giz_"+run_id+"_"+output_dir+anim_prams["suffix"]+"_%03d"%idump+"_"+plot_str+".dat"
                    np.savetxt(outfile,map)

    
    print("to mp4!")
    cmd = "ffmpeg -y -r 24 -hide_banner -loglevel quiet -stats -i "+pic_dir+"/sphplot"+run_id+output_dir+"%03d.png -c:v mpeg4 -q:v 1 "+movieDir+"/"+smooth_str+"_"+outp_plot+"_"+run_id+"_"+output_dir+anim_prams["suffix"]+".mp4"

    print(cmd)
    os.system(cmd)


if __name__ == '__main__':
    anim_prams = parse_input()
    animate(anim_prams)
