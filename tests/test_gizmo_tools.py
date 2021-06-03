from unittest import TestCase
import copy

from gtools import gizmo_tools
import numpy as np
import numpy.testing as npt

import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt


class Test_Functions(TestCase) :
    def test_load_interpolate_opacity(self) :
        self.assertIsInstance(gizmo_tools.load_interpolate_opacity(2.), np.ndarray)

        mu = np.arange(0.01, 100, 1000)
        opac = gizmo_tools.load_interpolate_opacity(mu)
        self.assertIsInstance(opac, np.ndarray)
        self.assertTrue(np.isfinite(opac) == len(opac))
        self.assertTrue(len(opac) == len(mu))
        self.assertTrue(opac.dtype == mu.dtype)

        f = gizmo_tools.dust_opacity_function
        self.assertIsInstance(f(2.), np.ndarray)
        opac2 = f(mu) * 0.000208908219  # unit conversion
        npt.assert_almost_equal(opac, opac2)

    def test_derive_opacities(self) :
        gizmo_tools.derive_opacities()
        lo = gizmo_tools.line_opacities
        self.assertIsInstance(lo, dict)
        self.assertSequenceEqual(list(lo.keys()), gizmo_tools.line_bases)
        for value in lo.values() :
            self.assertTrue(np.isfinite(value) == 1)

    def test_sort_nicely(self) :
        int_l = list(range(10000))
        int_l.reverse()
        str_l = list(map(str, int_l))
        gizmo_tools.sort_nicely(str_l)
        int_sorted = list(map(int, str_l))
        int_l.reverse()
        self.assertSequenceEqual(int_sorted, int_l)


class Test_Figs(TestCase) :
    def test_box_connected_two_axes(self) :
        fig, sp = plt.subplots(2, 1)
        oldpos0 = copy.copy(sp[0].get_position().get_points())
        oldpos1 = copy.copy(sp[1].get_position().get_points())
        gizmo_tools.box_connected_two_axes(sp[0], sp[1], [.5, 0.], [1., 1.5])
        gizmo_tools.box_connected_two_axes(sp[1], sp[0], [.5, -.5], [100., 1.], color='orange')
        npt.assert_equal(sp[0].get_position().get_points(), oldpos0)
        npt.assert_equal(sp[1].get_position().get_points(), oldpos1)


class DirTest(TestCase) :
    def test_get_dirs(self) :
        for f in [gizmo_tools.getDataDir,
                  gizmo_tools.getPicDir,
                  gizmo_tools.getTableDir,
                  gizmo_tools.getMovieDir]:
            self.assertIsInstance(f(),str)

    # how to test GetGizmoDir() ?


