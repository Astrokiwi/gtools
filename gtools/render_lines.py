# Import libraries to do our magic
import numpy as np
import sys
import os

from . import frame
from . import gizmo_tools
import os.path

import matplotlib.pyplot as plt


# with extinction
vranges = {"viewIRdust": [-5., 5.],
           "viewco1": [-20., 0.],
           "viewco2": [-20., 0.],
           "viewhcn1": [-20., 0.],
           "viewhcn2": [-20., 0.],
           #             "viewh2_1":[-20.,0.],
           #             "viewh2_2":[-20.,0.],
           #             "viewh2_3":[-20.,0.]
           "viewh2_1": [-10., 0.],
           "viewh2_3": [-10., 0.],
           "view12mic": [-20, 0.],
           "view8mic": [-20, 0.],
           "view850mic": [-20, 0.]
           }

# all lines
base_codes = ["co1", "co2", "hcn1", "hcn2", "h2_1", "h2_2", "h2_3", "IRdust", "12mic", "8mic", "850mic"]
viz_codes = ["view" + base_code for base_code in base_codes]
line_codes_all = [viz_code + "_000" for viz_code in viz_codes] + [viz_code + "_001" for viz_code in viz_codes]

# line_codes = ["viewh2_1_001","viewh2_3_001"]
#line_codes = ["view12mic_001", "view850mic_001", "viewco1_001"]
# line_codes = ["view12mic_001", "view850mic_001", "viewco1_001"]
line_codes = line_codes_all.copy()
rads = [100.] * 3

line_ids = [x[4:-4] for x in line_codes]
extinction_suffix = "extinction"
print(line_ids)

line_name_lookup = gizmo_tools.lineNamesFull.copy()
line_name_lookup["IRdust"] = r"$F_{IR}$"
line_labels = [line_name_lookup[line_id] for line_id in line_ids]

nlines = len(line_codes)  # len(rgb_image.labels)

cmap = "inferno"

scale = .7


def render_all(run_ids, output_dirs, idump):
    print("Rendering")
    snap_str = "{:03d}".format(idump)
    data_dir = gizmo_tools.getDataDir()
    pic_dir = gizmo_tools.getPicDir()

    itest = 0

    for run_id, output_dir in zip(run_ids, output_dirs) :
        gizmoDir = gizmo_tools.getGizmoDir(run_id)
        fullDir = os.path.join(gizmoDir,run_id,output_dir)
        infile = f"{fullDir}/snapshot_{idump:03d}.hdf5"
        print(infile)

        outfile = f"{pic_dir}/viztestv_{itest:03d}.png"

        # toplot = ["view12mic", "view850mic", "viewco1", "viewIRdust"]
        # data_ranges = [[-20., 0.], [-20., 0.], [-20., 0.], [-5, 5]]
        toplot = ["v12mic", "v850mic", "vco1", "v8mic"]
        data_ranges = [[-40.,40.]]*len(toplot)
        cmap = 'bwr'

        print(infile,outfile)
        m = frame.makesph_trhoz_frame_wrapped(infile, outfile,
                                    scale=.2,
                                    cmap=cmap,
                                    L=64,
                                    visibleAxes=True,
                                    pixsize=4,
                                    plot=toplot,
                                    rot=[0., 0.],
                                    views=['face'],
                                    data_ranges=data_ranges,
                                    return_maps=True,
                                    )

        print(len(m))

        outfile = f"{pic_dir}/viztestrotv_{itest:03d}.png"

        print(infile,outfile)
        m = frame.makesph_trhoz_frame_wrapped(infile, outfile,
                                    scale=.2,
                                    cmap=cmap,
                                    L=64,
                                    visibleAxes=True,
                                    pixsize=4,
                                    plot=toplot,
                                    rot=[0., np.pi/3],
                                    views=['face'],
                                    data_ranges=data_ranges,
                                    return_maps=True,
                                    )

        print(len(m))
        itest+=1

#
# def plot_lines_separately():
#     nx = 2
#     ny = nruns // nx
#     if ny * nx < nruns: ny += 1
#
#     for iline, line_code in enumerate(line_codes):
#         #     if line_code!="view": continue
#         line_label = rgb_image.line_label_dict[line_code]
#         print("plotting {}".format(line_label))
#         fig, sp = plt.subplots(ny, 1 + nx, figsize=(4. * nx + 0.25 * nx, 3. * ny),
#                                gridspec_kw={'width_ratios': [16] * nx + [1], 'height_ratios': [16] * ny}, squeeze=False)
#
#         for irow in range(ny - 1):
#             for icol in range(nx):
#                 sp[irow, icol].set_xticklabels([])
#             for icol in range(1, nx):
#                 sp[irow, icol].set_yticklabels([])
#             for icol in range(nx - 1):
#                 sp[irow, icol].get_shared_x_axes().join(sp[irow, icol], sp[irow + 1, icol])
#                 sp[irow, icol].get_shared_x_axes().join(sp[irow, icol], sp[irow, icol + 1])
#
#         for irow in range(ny):
#             for icol in range(nx):
#                 sp[irow, icol].set_aspect('equal', 'box')
#
#         plot_str = line_code
#         if line_codes.count(plot_str) > 0:
#             plot_str += "_%03d" % line_codes[:iline].count(line_code)
#
#         if line_code in vranges:
#             vrange = vranges[line_code]
#         else:
#             vrange = None
#
#         for irow, (run_id,output_dir) in enumerate(zip(run_ids,output_dirs)):
#             ix = irow % nx
#             iy = irow // nx
#             infile = "../data/smoothsum_giz_{}_{}_{:03d}_{}.dat".format(run_id, output_dir, idump, plot_str)
#             infile = [infile]
#
#             label = [line_label]
#             print("..subplot {}".format(infile[0]))
#
#             ax = sp[iy, ix]
#
#             mappable = rgb_image.produce_and_save_rgb_table_image(infile, ax=ax, rad=rads[iline], labels=label,
#                                                                   vrange=vrange, cmap='inferno')
#             ax.set_title(output_dir)
#             ax.set_xlabel("")
#             ax.set_ylabel("")
#             if ix == 0:
#                 plt.colorbar(mappable, cax=sp[iy, nx])
#
#         if nx * ny > nruns:
#             for iplot in range(nruns, nx * ny):
#                 ix = iplot % nx
#                 iy = iplot // nx
#                 sp[iy, ix].remove()
#
#         sp[ny - 1, 0].set_xlabel("pc")
#         sp[ny - 1, 0].set_ylabel("pc")
#
#         fig.tight_layout()
#         print(".finishing plot")
#
#         plt.savefig("../../figures/line_montage_{}.png".format(plot_str), dpi=150)
#
#
# def which_files_exist(run_ids, output_dirs, idump):
#     output_dirs_which_exist = []
#     for run_id,output_dir in zip(run_ids,output_dirs):
#         file_exists = True
#         for line_code in line_codes:
#             infile = "../data/smoothsum_giz_{}_{}_{:03d}_{}.dat".format(run_id, output_dir, idump, line_code)
#             if not os.path.isfile(infile):
#                 file_exists = False
#                 break
#         if file_exists:
#             output_dirs_which_exist.append(output_dir)
#     return output_dirs_which_exist
#
#
# def plot_lines_together(run_ids, output_dirs, idump, split=False):
#     ny = len(output_dirs)
#     nx = len(line_codes)
#     if split:
#         nx //= 2
#         ny *= 2
#     fw_inches = 3. * nx * scale
#     pixsize = 1
#     fig, sp = plt.subplots(ny + 1, nx, figsize=(fw_inches, scale * (3. * ny + 3. / 16.)),
#                            gridspec_kw={'width_ratios': [16] * nx, 'height_ratios': [16] * ny + [1]}, squeeze=False,
#                            constrained_layout=True)
#
#     for irow in range(ny):
#         for icol in range(nx):
#             sp[irow, icol].set_aspect('equal', 'box')
#
#     for iline, line_code in enumerate(line_codes):
#         for idir, (run_id,output_dir) in enumerate(zip(run_ids,output_dirs)):
#             infile = "../data/smoothsum_giz_{}_{}_{:03d}_{}.dat".format(run_id, output_dir, idump, line_code)
#             infile = [infile]
#
#             print("..subplot {}".format(infile[0]))
#             if split:
#                 ix = iline % (len(line_codes) // 2)
#                 iy = iline // (len(line_codes) // 2) + idir * 2
#             else:
#                 ix = iline
#                 iy = idir
#
#             ax = sp[iy, ix]
#             if line_code in vranges:
#                 vrange = vranges[line_code]
#             elif line_code[:-4] in vranges:
#                 vrange = vranges[line_code[:-4]]
#             else:
#                 vrange = None
#
#             mappable = rgb_image.produce_and_save_rgb_table_image(infile, ax=ax, rad=rads[iline], vrange=vrange,
#                                                                   cmap='inferno', printrange=True)  # ,gauss_rad=1.)
#
#             ax.set_title("")
#             if iy == ny - 1 or split and iy % 2 == 1:
#                 ax.set_xlabel("pc")
#             else:
#                 ax.set_xlabel("")
#                 ax.xaxis.set_visible(False)
#             if iy == 0:
#                 ax.set_title(line_labels[iline])
#                 plt.colorbar(mappable, cax=sp[ny, ix], orientation='horizontal')
#             if ix == 0:
#                 ax.set_ylabel(config3032.run_parameters[output_dir]["name"])
#             else:
#                 ax.yaxis.set_visible(False)
#                 ax.set_ylabel("")
#
#     print("Dumping png file")
#
#
#     L = 400
#     pos = sp[0, 0].get_position()
#     lpixx = (L / (pos.x1 - pos.x0))
#     my_dpi = int(np.floor(lpixx / fw_inches)) * pixsize
#
#     print(my_dpi)
#
#     plt.savefig("../../figures/lines_together_summary_{}_{}.png".format(extinction_suffix, idump), dpi=my_dpi)

if __name__ == '__main__':
    run_ids = ["rad_unary", "rad_prod", "rad_prod"]
    output_dirs = ["close", "rad_small_circ_earlier", "rad_small_ecc_earlier"]
    # run_ids = ["rad_prod"]
    # output_dirs = ["rad_small_ecc_earlier"]

    print("Running")
    idumps = [20]

    # if len(sys.argv) >= 2:
    for idump in idumps:
        render_all(run_ids,output_dirs,idump)
    # else:
    #     for idump in idumps:
    #         plot_lines_together(run_ids, which_files_exist(run_ids, output_dirs, idump), idump, split=False)
