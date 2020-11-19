import os


# Import libraries to do our magic
import tblib.pickling_support


import numpy as np
import pandas as pd
import pynbody
import h5py

import sys
import itertools
import json
import copy

from difflib import SequenceMatcher
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import colors
import matplotlib.pyplot as plt

from .sph_plotter import sph_plotter as sph_plotter_noisy
# from sph_plotter import sph_plotter
from . import output_suppressor
from . import preqs_config
from . import gizmo_tools




tblib.pickling_support.install()
sph_plotter = output_suppressor.suppressed(sph_plotter_noisy)
this_dir, this_filename = os.path.split(__file__)


class ExceptionWrapper(object):

    def __init__(self, ee):
        self.ee = ee
        __, __, self.tb = sys.exc_info()

    def re_raise(self):
        raise self.ee.with_traceback(self.tb)


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


dust_opacity_function = None


def load_interpolate_opacity(opac_mu,
                             opac_file=None):
    global dust_opacity_function
    if opac_file is None:
        opac_file = os.path.join(this_dir, "../prams/simple_dust_opacity.txt")
    if dust_opacity_function is None:
        dust_opacity_table = np.loadtxt(opac_file)
        dust_opacity_function = interpolate.interp1d(dust_opacity_table[:, 0], dust_opacity_table[:, 1])
    opacity = dust_opacity_function(opac_mu)
    opacity *= 0.000208908219  # convert to pc**2/solar mass for consistency
    return opacity


lines = ["co1", "co2", "hcn1", "hcn2", "h2_1", "h2_2", "h2_3", "12mic", "8mic", "850mic"]
line_wavelengths = [866.727, 433.438, 845.428, 422.796, 2.121, 28.18, 9.66, 12, 8, 850]
line_opacities = {line: load_interpolate_opacity(mu) for line, mu in zip(lines, line_wavelengths)}

# count "IRdust" as a line
lines = ["IRdust"] + lines
line_opacities.update({"IRdust": 0.1 * 0.000208908219})

densslice = 'densslice'
weightslice = 'weightslice'
zdensslice = 'zdensslice'
zweightslice = 'zweightslice'
vec2dslice_flat = 'vec2dslice_flat'
viewslice = 'viewslice'
minslice = 'minslice'
zminslice = 'zminslice'
vorinoislice = 'vorinoislice'
zvorinoislice = 'zvorinoislice'
maxslice = 'maxslice'
zmaxslice = 'zmaxslice'
maxdotslice = 'maxdotslice'
mindotslice = 'mindotslice'
sdevslice = 'sdevslice'
weightviewslice = 'weightviewslice'
zvec2dslice = 'zvec2dslice'
quiverslice_flat = 'quiverslice_flat'
zquiverslice = 'zquiverslice'

flat_choices = {'quantslice': [weightslice, zweightslice]
    , 'dslice': [densslice, zdensslice]
    , 'thisminslice': [minslice, zminslice]
    , 'vslice': [vorinoislice, zvorinoislice]
    , 'mslice': [maxslice, zmaxslice]
    , 'vec2dslice': [vec2dslice_flat, zvec2dslice]
    , 'quiverslice': [quiverslice_flat, zquiverslice]
                }

molecular_mass = 4. / (1. + 3. * .76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5. / 3. - 1.
boltzmann_cgs = 1.38066e-16

debug_mode = False

np.seterr(all='ignore')  # don't worry about bad logs etc - we want to propagate NaNs


def if_not_debug(f):
    def empty_function(*args, **kwargs):
        pass

    if debug_mode:
        return f
    else:
        return empty_function


@if_not_debug
def verboseprint(*args, **kwargs):
    print(*args, **kwargs)


# adapted from http://www.astrobetter.com/wiki/tiki-index.php?page=python_radial_profiles
def azimuthalNorm(image, center=None):
    """
    Normalize profile by max in each radial bin.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
    
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max() - x.min()) / 2.0, (x.max() - x.min()) / 2.0])

    r = np.hypot(x - center[0], y - center[1])

    rmax = np.max(r)
    nbins = int(np.ceil(rmax))

    ind_r = r.astype(int)

    bins_n = [np.max(image[ind_r == ii]) for ii in range(nbins)]

    for ii in range(nbins):
        image[ind_r == ii] = image[ind_r == ii] / bins_n[ii]

    return image


# class GadgetData:
#     pass

# wrapper for pcolormesh so it doesn't die if there's nothing finite to plot
def safe_pcolormesh(sp, *args, **kwargs):
    # args[2] is the actual map
    if np.sum(np.isfinite(args[2])) == 0:
        args = list(args)
        args[2][:, :] = kwargs["vmin"]
        args = tuple(args)
    return sp.pcolormesh(*args, shading='auto', **kwargs)


def makesph_plot(data, plane_keys,
                 m_p=None,
                 val_p=None,
                 plot_type='dslice',
                 h_p=None,
                 L=400,
                 width=1.,
                 corner='centered',
                 cblabel='VALUE',
                 vmin=0,
                 vmax=0,
                 cmap='iridis',
                 fig=None,
                 sp=None,
                 cbax=None,
                 z_p=None,
                 zslice=0.,
                 mask=None,
                 dolog=True,
                 planenorm=False,
                 circnorm=False,
                 plusminus=False,
                 diverging=False,
                 visibleAxes=True,
                 cbar2=None,
                 cmap2=None,
                 contour=False,
                 gaussian=None,
                 symLog=None,
                 cax_orientation='vertical',
                 centrecross=True,
                 flatPlot=True):
    x_p = data[plane_keys[0]]
    y_p = data[plane_keys[1]]
    if isinstance(val_p, str):
        val_p = data[val_p]

    cmap_label2 = cmap2
    cbax2 = cbar2

    this_cmap = copy.copy(plt.get_cmap(cmap))

    if cmap_label2:
        this_cmap2 = copy.copy(plt.get_cmap(cmap_label2))

    if not plusminus and not diverging:
        this_cmap.set_bad('black', 1.)
        this_cmap.set_under('black', 1.)

    n = x_p.size
    if mask is None:
        mask = np.ones(n)
    if corner == 'centered':
        corner = [-width / 2., -width / 2.]
    if m_p is None:
        m_p = np.ones(n)

    if plot_type in flat_choices:
        if flatPlot:
            mode = flat_choices[plot_type][0]
        else:
            mode = flat_choices[plot_type][1]
    else:
        mode = plot_type

    outmap = None

    if mode == weightslice:
        sph_map = sph_plotter.sph_weight(x_p, y_p, m_p, h_p, val_p, L, corner, width, mask, n)
    elif mode == densslice:
        sph_map = sph_plotter.sph_dense(x_p, y_p, m_p, h_p, L, corner, width, mask, n)
    elif mode == zdensslice:
        sph_map = sph_plotter.sph_dense_slice(x_p, y_p, m_p, h_p, L, corner, width, z_p, zslice, mask, n)
    elif mode == zweightslice:
        sph_map = sph_plotter.sph_weight_slice(x_p, y_p, m_p, h_p, val_p, L, corner, width, z_p, zslice, mask, n)
    elif mode == vec2dslice_flat or mode == quiverslice_flat:
        sph_map1 = sph_plotter.sph_weight(x_p, y_p, m_p, h_p, val_p[0], L, corner, width, mask, n)
        sph_map2 = sph_plotter.sph_weight(x_p, y_p, m_p, h_p, val_p[1], L, corner, width, mask, n)
    elif mode == zvec2dslice or mode == zquiverslice:
        sph_map1 = sph_plotter.sph_weight_slice(x_p, y_p, m_p, h_p, val_p[0], L, corner, width, z_p, zslice, mask, n)
        sph_map2 = sph_plotter.sph_weight_slice(x_p, y_p, m_p, h_p, val_p[1], L, corner, width, z_p, zslice, mask, n)
    elif mode == viewslice:
        zarg = np.argsort(z_p)
        sph_map = sph_plotter.sph_optical_depth_los_area(x_p, y_p, m_p, h_p, val_p[0], val_p[1], L, corner, width, z_p,
                                                         zarg, mask, n)

        total_luminosity = np.sum(sph_map) * (width * 2 * 3.086e+18 / L) ** 2 / 3.839e33  # check normalisation
        verboseprint(
            "log10(L/Lsun) = {}, mean flux = {} erg/s/cm**2, fullwidth = {} pc".format(np.log10(total_luminosity),
                                                                                       np.mean(sph_map), width * 2))
    elif mode == weightviewslice:
        zarg = np.argsort(z_p)
        #         threshold_flux = 1.e-5
        threshold_flux = 1.e-7
        print("HACKY: USING THRESHOLD FLUX OF ", threshold_flux)
        sph_map = sph_plotter.sph_optical_depth_los_weight_thresh(x_p, y_p, m_p, h_p, val_p[0], val_p[1], val_p[2], L,
                                                                  corner, width, z_p, zarg, mask, threshold_flux, n)


    elif mode == minslice:
        sph_map = sph_plotter.sph_min(x_p, y_p, h_p, val_p, L, corner, width, mask, n)
    elif mode == zminslice:
        sph_map = sph_plotter.sph_minslice(x_p, y_p, h_p, val_p, L, corner, width, mask, n)
    elif mode == vorinoislice:
        sph_map = sph_plotter.sph_vorinoi(x_p, y_p, m_p, h_p, val_p, L, corner, width, mask, n)
    elif mode == zvorinoislice:
        sph_map = sph_plotter.sph_vorinoi_slice(x_p, y_p, m_p, h_p, val_p, L, corner, width, z_p, zslice, mask, n)
    elif mode == maxslice:
        sph_map = sph_plotter.sph_max(x_p, y_p, m_p, h_p, val_p, L, corner, width, mask, n)
    elif mode == zmaxslice:
        sph_map = sph_plotter.sph_max_slice(x_p, y_p, m_p, h_p, val_p, L, corner, width, z_p, zslice, mask, n)
    elif mode == maxdotslice:
        sph_map = sph_plotter.sph_dot(x_p, y_p, val_p, L, corner, width, 0, mask, n)
    elif mode == mindotslice:
        sph_map = sph_plotter.sph_dot(x_p, y_p, val_p, L, corner, width, 1, mask, n)
    elif mode == sdevslice:
        sph_map = sph_plotter.sph_sdev(x_p, y_p, m_p, h_p, val_p, L, corner, width, mask, n)

    if mode == vec2dslice_flat or mode == zvec2dslice or mode == quiverslice_flat or mode == zquiverslice:
        norm_map = np.sqrt(sph_map1 ** 2 + sph_map2 ** 2)
        sph_map1 /= norm_map
        sph_map2 /= norm_map

        sph_map1 = sph_map1.T
        sph_map2 = sph_map2.T

        norm_map = norm_map.T

        step = width / L
        xmids = np.arange(corner[0] + step, corner[0] + width + step, step)
        ymids = np.arange(corner[1] + step, corner[1] + width + step, step)
        # sp.set_axis_bgcolor('black') # deprecated apparently?
        sp.set_facecolor('black')
        norm_map = np.log10(norm_map)
        xedges = np.arange(corner[0], corner[0] + width, width / L)
        yedges = np.arange(corner[1], corner[1] + width, width / L)

        if sp is not None:
            qv = safe_pcolormesh(sp, xedges, yedges, norm_map, cmap=this_cmap, vmin=vmin, vmax=vmax)
            if mode == quiverslice_flat or mode == zquiverslice:
                qv = sp.quiver(xmids, ymids, sph_map1, sph_map2, norm_map, headwidth=10., pivot='mid', cmap='jet')
            else:
                sp.streamplot(xmids, ymids, sph_map1, sph_map2)
            if cbax is not None:
                if visibleAxes:
                    cb = fig.colorbar(qv, label=cblabel, cax=cbax, orientation=cax_orientation)
                else:
                    cbax.set_axis_off()
        outmap = [norm_map, sph_map1, sph_map2]
    else:
        step = width / L
        xedges = np.arange(corner[0] + step / 2., corner[0] + width + step / 2., step)
        yedges = np.arange(corner[1] + step / 2., corner[1] + width + step / 2., step)
        if planenorm:
            sph_map[:, :] = (sph_map[:, :].T / np.max(sph_map[:, :], axis=1)).T
        if circnorm:
            sph_map = azimuthalNorm(sph_map)

        if plusminus:

            plusmap = sph_map.copy()
            minusmap = sph_map.copy()

            isplus = (sph_map > 0.)
            isminus = (sph_map < 0.)

            if np.sum(isplus) > 0:
                plusmap[isminus] = 0.
                if sp is not None:
                    mesh1 = safe_pcolormesh(sp, xedges, yedges, plusmap.T, cmap=this_cmap, vmin=vmin, vmax=vmax,
                                            norm=colors.LogNorm())
                    if cbax is not None:
                        if visibleAxes:
                            fig.colorbar(mesh1, label=cblabel, cax=cbax, orientation=cax_orientation)
                        else:
                            cbax.set_axis_off()

            if np.sum(isminus) > 0:
                minusmap[isplus] = 0.
                minusmap = -minusmap

                if sp is not None:
                    mesh2 = safe_pcolormesh(sp, xedges, yedges, minusmap.T, cmap=this_cmap2, vmin=vmin, vmax=vmax,
                                            norm=colors.LogNorm())
                    if cbax is not None:  # do want an error if cbax is given but not cbax2
                        if visibleAxes:
                            fig.colorbar(mesh2, label=cblabel, cax=cbax2, orientation=cax_orientation)
                        else:
                            cbax2.set_axis_off()
            outmap = [plusmap.T, minusmap.T]
        else:
            verboseprint(np.nanmin(sph_map), np.nanmax(sph_map))

            if gaussian is not None:
                if gaussian > 0.:
                    gaussian_pix = int(np.floor(L * gaussian / width))
                    if gaussian_pix >= 1:
                        sph_map = gaussian_filter(sph_map, gaussian)

            if dolog:
                sph_map = np.log10(sph_map)

            if sp is not None:
                if symLog is not None:
                    mesh = safe_pcolormesh(sp, xedges, yedges, sph_map.T, cmap=this_cmap, vmin=vmin, vmax=vmax,
                                           norm=colors.SymLogNorm(linthresh=symLog, linscale=symLog, vmin=vmin,
                                                                  vmax=vmax))
                else:
                    mesh = safe_pcolormesh(sp, xedges, yedges, sph_map.T, cmap=this_cmap, vmin=vmin, vmax=vmax)
            if cbax is not None:
                if visibleAxes:
                    cb = fig.colorbar(mesh, label=cblabel, cax=cbax, orientation=cax_orientation)
                else:
                    cbax.set_axis_off()

            if sp is not None:
                if contour:
                    sp.contour(xedges, yedges, gaussian_filter(sph_map.T, 5), vmin=vmin, vmax=vmax, colors='white',
                               levels=8, linewidths=.5)

            outmap = sph_map.T

    if sp is not None:
        sp.set_xlim([corner[0], corner[0] + width])
        sp.set_ylim([corner[1], corner[1] + width])
        if visibleAxes and centrecross:
            sp.plot([0], [0], '+g', markersize=10., markeredgewidth=1.)

        if data.binary_positions:

            verboseprint("binary at ", data.binary_positions)
            for ibinary, marker in enumerate(['x', '+']):
                bin_coords = [0, 0]
                for coord_id, plane_axis in enumerate(plane_keys):
                    if plane_axis == 'x':
                        bin_coords[coord_id] = data.binary_positions[ibinary][0]
                    elif plane_axis == 'y':
                        bin_coords[coord_id] = data.binary_positions[ibinary][1]
                    elif plane_axis == 'z':
                        bin_coords[coord_id] = data.binary_positions[ibinary][2]
                    elif plane_axis == 'r':
                        bin_coords[coord_id] = np.sqrt(
                            data.binary_positions[ibinary][1] ** 2 + data.binary_positions[ibinary][0] ** 2)
                sp.scatter(bin_coords[0], bin_coords[1], marker=marker)
        else:
            sp.plot([0], [0], '+g', markersize=10., markeredgewidth=1.)
    return outmap


def load_gadget(infile, plot_thing
                , centredens=False):
    f = h5py.File(infile, "r")

    header = f["/Header"]
    time = header.attrs.get("Time")
    time *= 0.9778e9  # to yr
    time /= 1.e6  # to Myr

    xyz = np.array(f["/PartType0/Coordinates"])  # kpc
    n = xyz.shape[0]
    data = pynbody.new(gas=n)
    data["xyz"] = xyz

    try:
        BH_data = f["/BH_binary"]
        Binary_pos_1 = BH_data.attrs.get("Binary_pos_1")
        Binary_pos_2 = BH_data.attrs.get("Binary_pos_2")
        if isinstance(Binary_pos_1, np.ndarray) & isinstance(Binary_pos_2, np.ndarray):
            data.binary_positions = [Binary_pos_1 * 1.e3, Binary_pos_2 * 1.e3]
        else:
            data.binary_positions = None
    #         verboseprint("Binary BH data loaded")
    except KeyError as e:
        data.binary_positions = None
    #         verboseprint("No Binary BH data found, skipping")

    data["m_p"] = np.array(f["/PartType0/Masses"])  # 10^10 msun
    data["m_p"] *= 1e+10  # 10^10 solar masses to solar masses
    data["h_p"] = np.array(f["/PartType0/SmoothingLength"])  # kpc

    need_to_load = set(plot_thing)
    if centredens:
        need_to_load.add("nH")
    need_to_load = preqs_config.process_preqs(need_to_load)

    if "id" in need_to_load:
        data["id_p"] = np.array(f["/PartType0/ParticleIDs"]).astype(int)

    if "age" in need_to_load:
        infile_split = infile.split("/")
        run_id = infile_split[-3]
        run_name = infile_split[-2]
        run_snapfile = infile_split[-1]
        run_isnap = int(run_snapfile[-8:-5])
        age_file = "../../data/age_{}_{}_{}.dat".format(run_id, run_name, run_isnap)
        data["age"] = time - np.loadtxt(age_file)

    if "pres" in need_to_load:
        data["pres"] = np.array(f["/PartType0/Pressure"])
    #         data.pres*=1.989e+33 # from internal units to dyne/cm*8*2

    if "arads" in need_to_load:
        data["arads"] = np.array(f["/PartType0/RadiativeAcceleration"])
        data["arads"] *= 3.24086617e-12  # to cm/s/s

    if "arad" in need_to_load:
        data["arad"] = np.sqrt(np.sum(data["arads"] ** 2, axis=1))

    if "accel" in need_to_load:
        data["accels"] = np.array(f["/PartType0/Acceleration"])
        data["accels"] *= 3.24086617e-12  # to cm/s/s
        data["accel"] = np.sqrt(np.sum(data["accels"] ** 2, axis=1))

    if "depth" in need_to_load:
        data["depth"] = np.array(f["/PartType0/AGNOpticalDepth"])  # Msun/kpc**2

    if "list" in need_to_load:
        data["list"] = np.arange(n)

    if "rand" in need_to_load:
        data["rand"] = np.random.random(n)

    if "nneigh" in need_to_load:
        data["nneigh"] = np.array(f["/PartType0/TrueNumberOfNeighbours"])

    if "temp" in need_to_load:
        data["u_p"] = np.array(f["/PartType0/InternalEnergy"])  # 1e10 erg/g
        data["u_p"] *= 1.e10  # to erg/g

        data["TK_p"] = (gamma_minus_one / boltzmann_cgs * (molecular_mass * proton_mass_cgs) * data["u_p"])

    if "vels" in need_to_load:
        data["vels"] = np.array(f["/PartType0/Velocities"])  # in km/s

    # doesn't work well
    if "dt" in need_to_load:
        data["dt_p"] = np.array(f["/PartType0/TimeStep"])
        data["dt_p"] *= 0.9778e9  # to yr

    # doesn't work well
    if "heat" in need_to_load:
        data["heat"] = np.array(f["/PartType0/AGNHeat"])
        data["heat"] *= 1e10 / 3.08568e+16  # to erg/s/g

    if "nH" in need_to_load:
        data["rho_p"] = np.array(f["/PartType0/Density"])
        data["rho_p"] *= 6.77e-22  # to g/cm**3
        data["nH_p"] = data["rho_p"] / (molecular_mass * proton_mass_cgs)

    if "tau" in need_to_load:
        data["tau"] = np.array(f["/PartType0/AGNDepth"])
    if "tau_2" in need_to_load:
        data["tau_2"] = np.array(f["/PartType0/AGNDepth_2"])
    if "tau_eff" in need_to_load:
        data["tau_eff"] = np.array(f["/PartType0/AGNDepth_Effective"])

    if "AGNI" in need_to_load:
        data["AGNI"] = np.array(f["/PartType0/AGNIntensity"])
        data["AGNI"] *= (1.989e53 / (
                0.9778e9 * 3.154e7) / 3.086e21 ** 2)  # convert from internal units (energy/Gyr-ish/kpc**2) to erg/s/cm**2

    if "AGNI2" in need_to_load:
        data["AGNI2"] = np.array(f["/PartType0/AGNIntensity_2"])
        data["AGNI2"] *= (1.989e53 / (
                0.9778e9 * 3.154e7) / 3.086e21 ** 2)  # convert from internal units (energy/Gyr-ish/kpc**2) to erg/s/cm**2

    if "AGNIeff" in need_to_load:
        data["AGNIeff"] = np.array(f["/PartType0/AGNIntensity_Effective"])
        data["AGNIeff"] *= (1.989e53 / (
                0.9778e9 * 3.154e7) / 3.086e21 ** 2)  # convert from internal units (energy/Gyr-ish/kpc**2) to erg/s/cm**2

    if "table" in need_to_load:
        verboseprint("Load dust tables")

        #         table_date="060319" # used in paper - not all intensities are there
        table_date = "130720"
        table_res = "0.1"
        coolheat_dir = os.path.join(this_dir, "../coolheat_tab_marta")
        cloudy_table = gizmo_tools.cloudy_table(table_date, table_res, coolheat_dir)
        data["flux_p"] = np.array(f["/PartType0/AGNIntensity"])  # energy per surface area per time
        data["flux_p"] *= 1.989e+53 / 3.086e21 ** 2 / 3.08568e+16

        table_particles = pd.DataFrame()
        table_particles["nH"] = data["nH_p"]
        table_particles["temp"] = data["TK_p"]
        table_particles["AGNIntensity"] = data["flux_p"]
        table_particles["AGNDepth"] = data["tau"]
        # TODO: correct optical depths etc for binary
        #         if "tau_eff" in data:
        #             table_particles["AGNDepth"] = data["tau_eff"]
        #         else:
        #             table_particles["AGNDepth"] = data["tau"]

        verboseprint("Calculating dust/cooling/heating properties from table")
        cloudy_table.interp(table_particles)

    if "col" in need_to_load:
        if "/PartType0/AGNColDens" in f:
            data["coldens"] = np.array(f["/PartType0/AGNColDens"])  # Msun/kpc**2
            data["coldens"] *= (1.989e+43 / 3.086e+21 ** 2)  # to g/cm**2
            data["coldens"] /= (molecular_mass * proton_mass_cgs)  # N in cm**(-2)
        else:
            data["coldens"] = table_particles["column_out"]

    if "tdust" in need_to_load:
        data["dustTemp"] = table_particles["dustT"]

    for line in lines[1:]:
        if line in need_to_load:
            data[line] = table_particles["line_" + line]
            if line in ["12mic", "8mic", "850mic"]:  # these are given in erg/cm**3/s, need to convert to erg/g/s
                data[line] /= data["nH_p"]
        if line + "m" in need_to_load:
            data[line + "m"] = table_particles["line_" + line] * data["m_p"] * 1.9891e33 / 9.52140614e36 / (
                    4. * np.pi)  # erg/s/g to erg/s, extra factor for pc**2 to ster cm**2, output is erg/s/cm**2/ster

    if "dg" in need_to_load:
        data["dg"] = table_particles["dg"]

    if "dust" in need_to_load:
        data["dust"] = data["dg"] * data["m_p"]

    if "emit" in need_to_load:
        data["emissivity"] = 5.67e-5 * data["dustTemp"] ** 4. * data["dg"] / np.nanmax(data["dg"])  # erg/s/cm^2

    if "opac" in need_to_load:
        data["opac"] = np.array(f["/PartType0/AGNOpacity"])  # internal units: kpc**2/1e10 solar mass
        data["opac"] *= 0.478679108  # to cm**2/g

    if ("view" in need_to_load or "vlos" in need_to_load or any(
            x in need_to_load for x in ["view" + line for line in lines])
            or "dusttau" in need_to_load):
        if "view" in need_to_load:
            data["brightness"] = 5.67e-5 * data["dustTemp"] ** 4. * data["dg"] / np.nanmax(data["dg"])  # erg/s/cm^2

        verboseprint("faking opacity")
        opacity = 65.2  # cm^2/g, somewhat arbitrary
        verboseprint("Broad IR opacity is now ", opacity, " cm^2/g")

        opacity *= 0.000208908219  # convert to pc**2/solar mass for consistency
        data["opac"] = np.full(n, opacity)
        data["opac"] *= data["dg"] / np.nanmax(data["dg"])  # take into account dust fraction

        if "dusttau" in need_to_load:
            data["dusttau"] = data["opac"] * data["m_p"]
        for line in lines[1:]:
            if not "view" + line in need_to_load:
                continue
            data[line + "opac"] = np.full(n, line_opacities[line])
            data[line + "brightness"] = data[line + "m"] / data[
                line + "opac"]  # erg/s/cm^2 - multiply by opacity to get actual emission
        #             data.__dict__[line+"brightness"]=data.__dict__[line+"m"] # erg/s, gets SPH smoothed to get erg/s/cm**2

        if "view" in need_to_load:
            # sputtered = (data.dustTemp>2.5e3) # or something, super arbitrary
            sputtered = (data["dustTemp"] > 1.e5)  # or something, super arbitrary
            data["brightness"][sputtered] = 0.
            data["opac"][sputtered] = 0.
    if any(x in need_to_load for x in ("IRdustm", "IRdust", "IRdustopac", "IRdustbrightness", "viewIRdust")):
        #         data["dustTemp"] = table_particles["dustT"]
        data["IRdustbrightness"] = 5.67e-5 * data["dustTemp"] ** 4. * data["dg"] / np.nanmax(data["dg"])
        data["IRdustopac"] = np.full(n, line_opacities["IRdust"])
        data["IRdustm"] = data["IRdustbrightness"] * data["IRdustopac"]
        data["IRdust"] = data["IRdustm"] / (data["m_p"] * 1.9891e33 / 9.52140614e36 / (4. * np.pi))

    if "rad0" in need_to_load:
        data["rad0"] = np.load("rad0.npy")
        data["rad0"] = data["rad0"][data["id_p"] - 1]

    return time, data


def pack_dicts():
    with open(os.path.join(this_dir, "plot_defaults.json"), 'r') as f:
        plot_config = json.load(f)

    # assign defaults
    for d in plot_config.values():
        if isinstance(d["range"], str):
            d["range"] = plot_config[d["range"]]["range"]
    return plot_config


def load_process_gadget_data(infile, rot, plot_thing, plot_config, centredens=False, ringPlot=False, flatPlot=False,
                             maskbounds=None):
    # ,opac_mu=None):
    need_to_load = list(plot_thing)
    if maskbounds:
        need_to_load.append(maskbounds[0])

    time, data = load_gadget(infile, need_to_load, centredens=centredens)  # ,opac_mu=opac_mu)

    n = len(data)

    # corotate with binaries?
    #     if data.binary_positions is not None:
    #       x_bin = data.binary_positions[0][0]-data.binary_positions[1][0]
    #       y_bin = data.binary_positions[0][1]-data.binary_positions[1][1]
    #       rot[0] = rot[0] - math.atan2(y_bin,x_bin)

    # convert to pc
    data["pos"] = data["xyz"] * 1.e3
    data["h_p"] *= 1.e3

    x = data["x"]
    y = data["y"]
    z = data["z"]

    # rotate
    if rot[0] != 0. or rot[1] != 0.:
        xr = x * np.cos(rot[0]) - y * np.sin(rot[0])
        yr = x * np.sin(rot[0]) + y * np.cos(rot[0])
        x = xr

        y = yr * np.cos(rot[1]) - z * np.sin(rot[1])
        z = yr * np.sin(rot[1]) + z * np.cos(rot[1])

        data["x"] = x
        data["y"] = y
        data["z"] = z

        if data.binary_positions:
            data.binary_positions_rot = [None, None]
            for ibin, bpos in enumerate(data.binary_positions):
                xrb = bpos[0] * np.cos(rot[0]) - bpos[1] * np.sin(rot[0])
                yrb = bpos[0] * np.sin(rot[0]) + bpos[1] * np.cos(rot[0])

                xb = xrb
                yb = yrb * np.cos(rot[1]) - bpos[2] * np.sin(rot[1])
                zb = yrb * np.sin(rot[1]) + bpos[2] * np.cos(rot[1])

                data.binary_positions_rot[ibin] = np.array([xb, yb, zb]) * 1.e3
    else:
        if data.binary_positions:
            data.binary_positions_rot = [data.binary_positions[0] * 1.e3, data.binary_positions[1] * 1.e3]

    if (any(x in plot_thing for x in
            ["vels", "vmag", "vel_2d", "vel_x", "vel_y", "vel_z", "vel_r", "vel_a", "vthin", "vlos"] + ["v" + line for
                                                                                                        line in
                                                                                                        lines] + [
                "dv" + line for line in lines] + ["vels" + line for line in
                                                  lines])):  # + ["view"+line for line in lines]
        # vel_mag = np.sqrt(np.sum(data.vels[:,:]**2,1))
        data["vel_a"] = (-x * data["vels"][:, 1] + y * data["vels"][:, 0]) / np.sqrt(x ** 2 + y ** 2)
        data["velr"] = (x * data["vels"][:, 0] + y * data["vels"][:, 1] + z * data["vels"][:, 2]) / np.sqrt(
            x ** 2 + y ** 2 + z ** 2)
        data["vel_x"] = data["vels"][:, 0]
        data["vel_y"] = data["vels"][:, 1]
        data["vel_z"] = data["vels"][:, 2]
        data["vmag"] = np.sqrt(np.sum(data["vels"] ** 2, axis=1))
        if rot[0] != 0. or rot[1] != 0.:
            if ringPlot:
                raise NotImplementedError()
            else:
                vyr = data["vel_y"] * np.cos(rot[1]) - data["vel_z"] * np.sin(rot[1])
                vzr = data["vel_y"] * np.sin(rot[1]) + data["vel_z"] * np.cos(rot[1])
                data["vel_z"] = vzr

                data["vel_y"] = data["vel_x"] * np.sin(rot[0]) + vyr * np.cos(rot[0])
                data["vel_x"] = data["vel_x"] * np.cos(rot[0]) - vyr * np.sin(rot[0])

                data["vel2d"] = data["vel_x"]
        else:
            if ringPlot:
                data["vel2d"] = (x * data["vel_x"] + y * data["vel_y"]) / np.sqrt(x ** 2 + y ** 2)
            else:
                data["vel2d"] = data["vel_x"]

    if any(x in plot_thing for x in ["arads", "arad_x", "arad_y", "arad_z", "arad_2d"]):
        data["arad_x"] = data["arads"][:, 0]
        data["arad_y"] = data["arads"][:, 1]
        data["arad_z"] = data["arads"][:, 2]
        if rot[0] != 0. or rot[1] != 0.:
            raise NotImplementedError()  # can't be arsed
        else:
            if ringPlot:
                data["arad2d"] = (x * data["arad_x"] + y * data["arad_y"]) / np.sqrt(x ** 2 + y ** 2)
            else:
                data["arad2d"] = data["arad_x"]

    if any(x in plot_thing for x in ["vthin", "vlos"] + ["v" + line for line in lines] + ["dv" + line for line in
                                                                                          lines]):  # + ["view"+line for line in lines]
        if ringPlot:
            raise Exception("can't do optically thin line of sight velocity while integrating around ring")

        # rotation already done above
        data["vthin"] = data["vel_z"]
        data["vythin"] = data["vel_y"]

    # set x coordinate - cartesian or cylindrical ("ringplot")
    if ringPlot:
        if not flatPlot:
            raise Exception("ring==true requires flat==true")
        data["r"] = np.sqrt(x ** 2 + y ** 2)
    else:
        data["r"] = x

    deep_face = z
    deep_side = y

    if maskbounds:
        if maskbounds[0] in plot_config:
            v = data[plot_config[maskbounds[0]]]
            mask = (v > maskbounds[1]) & (v < maskbounds[2])
        else:
            raise Exception("mask value {} not found in plot_config".format(maskbounds[0]))
    else:
        # mask out non-gas - currently everything is gas
        mask = np.full(n, True, dtype=bool)

    # flat weighting - dummy value required for some functions because I'm not qwarging properly yet
    n_ones = np.ones(n)

    return time, data, x, y, z, deep_face, deep_side, mask, n, n_ones


def makesph_trhoz_frame(*args, **kwargs):
    try:
        return makesph_trhoz_frame_wrapped(*args, **kwargs)
    except Exception as e:
        print("returning exception")
        wropped = ExceptionWrapper(e)
        return wropped


def makesph_trhoz_frame_wrapped(infile, outfile,
                                scale=.7,
                                cmap="viridis",
                                L=256,
                                ring=False,
                                flat=False,
                                planenorm=False,
                                visibleAxes=True,
                                subsample=1,
                                pixsize=None,
                                plot=['dens', 'temp'],
                                rot=[0., 0.],
                                views=['face', 'side'],
                                centredens=False,
                                centrecom=False,
                                vorinoi=False,
                                maxmode=False,
                                dotmode=None,
                                titlesuffix="",
                                maskbounds=None,
                                data_ranges=None,
                                gaussian=None,
                                return_maps=False
                                # ,opac_mu=None
                                ):
    if sys.version[0] == '2':
        raise Exception("Requires Python 3.X. Current version:" + sys.version)

    plot_thing = plot
    flatPlot = flat
    ringPlot = ring

    if return_maps:
        out_maps = []

    if dotmode != 'max' and dotmode != 'min':
        dotmode = None

    if L % subsample != 0:
        raise Exception("subsample might divide evenly into L")

    if not pixsize:
        pixsize = subsample

    if len(rot) != 2:
        raise Exception("rot needs to be [theta,phi]")

    if not all(view == 'face' or view == 'side' for view in views):
        raise Exception("Views must be an array of size one or two containing only 'face' or 'side'")
    if not 'side' in views:
        ringPlot = False  # don't care about ring plot if there's no side plot anyway

    cols = len(views)
    if cols != 1 and cols != 2:
        raise Exception("len(views)=1 or =2")

    if data_ranges is not None:
        if len(data_ranges) != len(plot):
            raise Exception("Need to give data ranges for ALL plots or none")

    plot_config = pack_dicts()

    nrows = len(plot_thing)

    if visibleAxes:
        fig, ax = plt.subplots(nrows, 2 * cols, gridspec_kw={'width_ratios': ([1, 16, 16, 1])[0:2 * cols]}, squeeze=False)
        cbaxleft_index = 0
        spleft_index = 1
        spright_index = 2
        cbaxright_index = 3
    else:
        fig, ax = plt.subplots(nrows, cols, squeeze=False)
        cbaxleft_index = 0
        spleft_index = 0
        spright_index = 1
        cbaxright_index = 1

    if not isinstance(infile, list):
        verboseprint("Loading", infile)
        time, data, x, y, z, deep_face, deep_side, mask, n, n_ones = load_process_gadget_data(infile, rot, plot_thing,
                                                                                              plot_config, centredens,
                                                                                              ringPlot, flatPlot,
                                                                                              maskbounds)  # ,opac_mu=opac_mu)
        #         time,data = load_gadget(infile,need_to_load,centredens=centredens)
        verboseprint("Plotting", infile, ", t=%.4f Myr" % time)

        if visibleAxes:
            if time < 1.e-3:
                fig.suptitle(r"$T=" + ("%.4f" % (time * 1.e3)) + "$ kyr" + titlesuffix)
            else:
                fig.suptitle(r"$T=" + ("%.4f" % time) + "$ Myr" + titlesuffix)
    #         fig.suptitle(infile)

    if visibleAxes:
        fw_inches = 5. * cols
    else:
        fw_inches = 4. * cols
    fig.set_figwidth(fw_inches)
    fig.set_figheight(4. * nrows)

    # nerrors=0

    # do all subplots, calculating the full SPH smoothing each time
    for irow in range(nrows):

        if isinstance(infile, list):
            thisfile = infile[irow]
            verboseprint("Loading", thisfile)
            time, data, x, y, z, deep_face, deep_side, mask, n, n_ones = load_process_gadget_data(infile, rot,
                                                                                                  plot_thing,
                                                                                                  plot_config,
                                                                                                  centredens,
                                                                                                  ringPlot, flatPlot,
                                                                                                  maskbounds)

            verboseprint("Plotting", thisfile, ", t=%.4f Myr" % time)

        # Process the tables we just made
        if not plot_thing[irow] in plot_config:
            errstr = "{} is not a valid plot type\n".format(plot_thing[irow])
            labelstrs = list(plot_config.keys())
            matchRatios = list(map(similar, itertools.repeat(plot_thing[irow]), labelstrs))
            bestFit = np.argmax(matchRatios)
            errstr += "Valid plot types:\n"
            for i in range(len(labelstrs)):
                if i != bestFit:
                    if matchRatios[i] > .5:
                        labelstrs[i] += " - possible fit?"
            labelstrs[bestFit] += " - best match"
            errstr += "\n".join(labelstrs)
            raise Exception(errstr)

        thisPlotLabel = plot_config[plot_thing[irow]]['label']
        if isinstance(thisPlotLabel, list):
            if flatPlot:
                thisPlotLabel = thisPlotLabel[0]
            else:
                thisPlotLabel = thisPlotLabel[1]

        if data_ranges is not None:
            thisPlotRanges = data_ranges[irow] * 2
        else:
            thisPlotRanges = plot_config[plot_thing[irow]]['range']

        if isinstance(gaussian, list):
            plot_gaussian = gaussian[irow]
        else:
            plot_gaussian = gaussian

        # physical coordinates of region to plot, in pc
        if isinstance(scale, list):
            width = scale[irow] * 2.
        else:
            width = scale * 2.

        if centredens or centrecom:
            if centredens and centrecom:
                raise Exception("Can't set both centredens and centrecom")
            if centredens:
                i_maxdens = np.argmax(data["nH_p"])
                corners = [x[i_maxdens] - width / 2, y[i_maxdens] - width / 2.]
                corners_side = [data["r"][i_maxdens] - width / 2., z[i_maxdens] - width / 2.]
            if centrecom:
                x_com = np.mean(x)
                y_com = np.mean(y)
                corners = [x_com - width / 2., y_com - width / 2.]
                r2d_com = np.mean(data["r"])
                z_com = np.mean(z)
                corners_side = [r2d_com - width / 2., z_com - width / 2.]
        else:
            # centre on 0,0
            corners = [-width / 2., -width / 2.]
            if ringPlot:
                if not flatPlot:
                    raise Exception("ring==true requires flat==true")
                corners_side = [0., -width / 2.]
            else:
                corners_side = corners

        if dotmode == 'max':
            thisSliceType = maxdotslice
        elif dotmode == 'min':
            thisSliceType = mindotslice
        else:
            thisSliceType = plot_config[plot_thing[irow]]['slice']

        thisDoLog = plot_config[plot_thing[irow]]['log']
        thisDiverging = plot_config[plot_thing[irow]]['diverging']
        thisSymLog = plot_config[plot_thing[irow]]['symlog']
        if not thisSymLog:
            thisSymLog = None

        thisMass = data["m_p"]

        plusminus = False  # NEVER USED

        if plot_thing[irow] in plot_config:
            if 'data' not in plot_config[plot_thing[irow]]:
                thisPlotQuantityFace = None
                thisPlotQuantitySide = None
            else:
                plotCommand = plot_config[plot_thing[irow]]['data']
                if type(plotCommand) is str:
                    thisPlotQuantityFace = data[plotCommand]
                    thisPlotQuantitySide = thisPlotQuantityFace
                elif type(plotCommand) is list:
                    if thisSliceType == viewslice or thisSliceType == 'vec2dslice' or thisSliceType == 'quiverslice':
                        if len(plotCommand) != 4:
                            raise Exception("{} must have length 4 for view or vec2d slice".format(plotCommand))
                        thisPlotQuantityFace = [data[plotCommand[0]], data[plotCommand[1]]]
                        thisPlotQuantitySide = [data[plotCommand[2]], data[plotCommand[3]]]
                    elif thisSliceType == weightviewslice:
                        if len(plotCommand) != 6:
                            raise Exception("{} must have length 6 for weighted view slice".format(plotCommand))
                        thisPlotQuantityFace = [data[plotCommand[0]], data[plotCommand[1]], data[plotCommand[2]]]
                        thisPlotQuantitySide = [data[plotCommand[3]], data[plotCommand[4]], data[plotCommand[5]]]
                    else:
                        raise Exception("{} can't be an array except for view and vec2d slices".format(plotCommand))
                else:
                    raise Exception("{} is not a valid plot format".format(plotCommand))

        cbar2_axes = [None, None]

        this_cmap = cmap

        this_cmap2 = None

        row_axes = [ax[irow, spleft_index], ax[irow, cbaxleft_index]]
        if cols == 2:
            row_axes += ax[irow, spright_index], ax[irow, cbaxright_index]

        # actually do the plot
        for icol, view in enumerate(views):
            if view == 'face':

                outmap = makesph_plot(data, 'xy'
                                      , m_p=thisMass
                                      , val_p=thisPlotQuantityFace
                                      , plot_type=thisSliceType
                                      , h_p=data["h_p"]
                                      , L=L
                                      , width=width
                                      , corner=corners
                                      , cblabel=thisPlotLabel
                                      , vmin=thisPlotRanges[0]
                                      , vmax=thisPlotRanges[1]
                                      , cmap=this_cmap
                                      , fig=fig
                                      , sp=row_axes[icol * 2]
                                      , cbax=row_axes[icol * 2 + 1]
                                      , z_p=deep_face
                                      , zslice=0.
                                      , mask=mask
                                      , dolog=thisDoLog
                                      , cmap2=this_cmap2
                                      , circnorm=planenorm
                                      , cbar2=cbar2_axes[icol]
                                      , plusminus=plusminus
                                      , visibleAxes=visibleAxes
                                      , diverging=thisDiverging
                                      , gaussian=plot_gaussian
                                      , symLog=thisSymLog)

            #                 if data.binary_positions:
            #                     row_axes[icol*2].scatter([data.binary_positions_rot[0][0]],[data.binary_positions_rot[0][1]],marker='x')
            #                     row_axes[icol*2].scatter([data.binary_positions_rot[1][0]],[data.binary_positions_rot[1][1]],marker='+')
            elif view == 'side':
                if ringPlot:
                    sideview = 'rz'
                else:
                    sideview = 'yz'
                outmap = makesph_plot(data, sideview
                                      , vmin=thisPlotRanges[2]
                                      , vmax=thisPlotRanges[3]
                                      , val_p=thisPlotQuantitySide
                                      , corner=corners_side
                                      , m_p=thisMass
                                      , plot_type=thisSliceType
                                      , h_p=data["h_p"]
                                      , L=L
                                      , width=width
                                      , cblabel=thisPlotLabel
                                      , cmap=this_cmap
                                      , fig=fig
                                      , sp=row_axes[icol * 2]
                                      , cbax=row_axes[icol * 2 + 1]
                                      , z_p=deep_side
                                      , zslice=0.
                                      , mask=mask
                                      , dolog=thisDoLog
                                      , cmap2=this_cmap2
                                      , circnorm=planenorm
                                      , cbar2=cbar2_axes[icol]
                                      , plusminus=plusminus
                                      , visibleAxes=visibleAxes
                                      , diverging=thisDiverging
                                      , gaussian=plot_gaussian
                                      , symLog=thisSymLog)
            if return_maps:
                out_maps.append(outmap)

    #         if ( "heat" in plot_thing and plot_thing[irow]!='heat' ):
    #             for visax in [ax[irow,cbax2left_index],ax[irow,cbax2right_index]]:
    #                 visax.set_frame_on(False)
    #                 visax.axes.get_yaxis().set_visible(False)
    #                 visax.axes.get_xaxis().set_visible(False)

    if visibleAxes:
        if cols == 2:
            for iax in range(nrows):
                this_ax = ax[iax, 1]
                this_ax.yaxis.tick_right()
                this_ax.yaxis.set_visible(False)
            for iax in range(nrows):
                this_ax = ax[iax, 0]
                this_ax.yaxis.tick_left()
                this_ax.yaxis.set_label_position("left")

        else:
            for iax in range(nrows):
                this_ax = ax[iax, 0]
                this_ax.yaxis.tick_left()
                this_ax.yaxis.set_label_position("left")
            for iax in range(nrows):
                this_ax = ax[iax, 1]
                this_ax.yaxis.tick_right()
                this_ax.yaxis.set_label_position("right")

    else:
        plt.axis('off')

    #     print(infile,nerrors)

    if visibleAxes:
        if nrows == 2:
            fig.subplots_adjust(left=0.07, hspace=.07, bottom=.05, top=.95)
        else:
            #             fig.subplots_adjust(left=0.07,hspace=.07,bottom=.05,top=.9)
            fig.subplots_adjust(left=0.12, hspace=.12, bottom=.1, top=.9)
    else:
        fig.subplots_adjust(left=0.0, hspace=.0, top=1., bottom=.0, right=1., wspace=0.)

    pos = ax[0, 1].get_position()
    lpixx = (L / (pos.x1 - pos.x0))
    my_dpi = int(np.floor(lpixx / fw_inches)) * pixsize

    fig.savefig(outfile, dpi=my_dpi)
    plt.close()
    if return_maps:
        return out_maps
    else:
        return None
