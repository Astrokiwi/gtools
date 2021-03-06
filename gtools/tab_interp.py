# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _tab_interp
else:
    import _tab_interp

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


class coolHeatDust(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    dCool = property(_tab_interp.coolHeatDust_dCool_get, _tab_interp.coolHeatDust_dCool_set)
    dHeat = property(_tab_interp.coolHeatDust_dHeat_get, _tab_interp.coolHeatDust_dHeat_set)
    dustT = property(_tab_interp.coolHeatDust_dustT_get, _tab_interp.coolHeatDust_dustT_set)
    opac_abs = property(_tab_interp.coolHeatDust_opac_abs_get, _tab_interp.coolHeatDust_opac_abs_set)
    opac_scat = property(_tab_interp.coolHeatDust_opac_scat_get, _tab_interp.coolHeatDust_opac_scat_set)
    dg = property(_tab_interp.coolHeatDust_dg_get, _tab_interp.coolHeatDust_dg_set)
    column_out = property(_tab_interp.coolHeatDust_column_out_get, _tab_interp.coolHeatDust_column_out_set)
    arad = property(_tab_interp.coolHeatDust_arad_get, _tab_interp.coolHeatDust_arad_set)
    line_co1 = property(_tab_interp.coolHeatDust_line_co1_get, _tab_interp.coolHeatDust_line_co1_set)
    line_co2 = property(_tab_interp.coolHeatDust_line_co2_get, _tab_interp.coolHeatDust_line_co2_set)
    line_hcn1 = property(_tab_interp.coolHeatDust_line_hcn1_get, _tab_interp.coolHeatDust_line_hcn1_set)
    line_hcn2 = property(_tab_interp.coolHeatDust_line_hcn2_get, _tab_interp.coolHeatDust_line_hcn2_set)
    line_h2_1 = property(_tab_interp.coolHeatDust_line_h2_1_get, _tab_interp.coolHeatDust_line_h2_1_set)
    line_h2_2 = property(_tab_interp.coolHeatDust_line_h2_2_get, _tab_interp.coolHeatDust_line_h2_2_set)
    line_h2_3 = property(_tab_interp.coolHeatDust_line_h2_3_get, _tab_interp.coolHeatDust_line_h2_3_set)
    line_12m = property(_tab_interp.coolHeatDust_line_12m_get, _tab_interp.coolHeatDust_line_12m_set)
    line_8m = property(_tab_interp.coolHeatDust_line_8m_get, _tab_interp.coolHeatDust_line_8m_set)
    line_850m = property(_tab_interp.coolHeatDust_line_850m_get, _tab_interp.coolHeatDust_line_850m_set)

    def __init__(self):
        _tab_interp.coolHeatDust_swiginit(self, _tab_interp.new_coolHeatDust())
    __swig_destroy__ = _tab_interp.delete_coolHeatDust

# Register coolHeatDust in _tab_interp:
_tab_interp.coolHeatDust_swigregister(coolHeatDust)

class coolHeatDustArray(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    dCool = property(_tab_interp.coolHeatDustArray_dCool_get, _tab_interp.coolHeatDustArray_dCool_set)
    dHeat = property(_tab_interp.coolHeatDustArray_dHeat_get, _tab_interp.coolHeatDustArray_dHeat_set)
    dustT = property(_tab_interp.coolHeatDustArray_dustT_get, _tab_interp.coolHeatDustArray_dustT_set)
    opac_abs = property(_tab_interp.coolHeatDustArray_opac_abs_get, _tab_interp.coolHeatDustArray_opac_abs_set)
    opac_scat = property(_tab_interp.coolHeatDustArray_opac_scat_get, _tab_interp.coolHeatDustArray_opac_scat_set)
    dg = property(_tab_interp.coolHeatDustArray_dg_get, _tab_interp.coolHeatDustArray_dg_set)
    column_out = property(_tab_interp.coolHeatDustArray_column_out_get, _tab_interp.coolHeatDustArray_column_out_set)
    arad = property(_tab_interp.coolHeatDustArray_arad_get, _tab_interp.coolHeatDustArray_arad_set)
    line_co1 = property(_tab_interp.coolHeatDustArray_line_co1_get, _tab_interp.coolHeatDustArray_line_co1_set)
    line_co2 = property(_tab_interp.coolHeatDustArray_line_co2_get, _tab_interp.coolHeatDustArray_line_co2_set)
    line_hcn1 = property(_tab_interp.coolHeatDustArray_line_hcn1_get, _tab_interp.coolHeatDustArray_line_hcn1_set)
    line_hcn2 = property(_tab_interp.coolHeatDustArray_line_hcn2_get, _tab_interp.coolHeatDustArray_line_hcn2_set)
    line_h2_1 = property(_tab_interp.coolHeatDustArray_line_h2_1_get, _tab_interp.coolHeatDustArray_line_h2_1_set)
    line_h2_2 = property(_tab_interp.coolHeatDustArray_line_h2_2_get, _tab_interp.coolHeatDustArray_line_h2_2_set)
    line_h2_3 = property(_tab_interp.coolHeatDustArray_line_h2_3_get, _tab_interp.coolHeatDustArray_line_h2_3_set)
    line_12mic = property(_tab_interp.coolHeatDustArray_line_12mic_get, _tab_interp.coolHeatDustArray_line_12mic_set)
    line_8mic = property(_tab_interp.coolHeatDustArray_line_8mic_get, _tab_interp.coolHeatDustArray_line_8mic_set)
    line_850mic = property(_tab_interp.coolHeatDustArray_line_850mic_get, _tab_interp.coolHeatDustArray_line_850mic_set)

    def __init__(self):
        _tab_interp.coolHeatDustArray_swiginit(self, _tab_interp.new_coolHeatDustArray())
    __swig_destroy__ = _tab_interp.delete_coolHeatDustArray

# Register coolHeatDustArray in _tab_interp:
_tab_interp.coolHeatDustArray_swigregister(coolHeatDustArray)

class AGN_heat_table(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    agn_heat_tab = property(_tab_interp.AGN_heat_table_agn_heat_tab_get, _tab_interp.AGN_heat_table_agn_heat_tab_set)
    agn_cool_tab = property(_tab_interp.AGN_heat_table_agn_cool_tab_get, _tab_interp.AGN_heat_table_agn_cool_tab_set)
    agn_dust_tab = property(_tab_interp.AGN_heat_table_agn_dust_tab_get, _tab_interp.AGN_heat_table_agn_dust_tab_set)
    agn_dg_tab = property(_tab_interp.AGN_heat_table_agn_dg_tab_get, _tab_interp.AGN_heat_table_agn_dg_tab_set)
    agn_opac_scat_tab = property(_tab_interp.AGN_heat_table_agn_opac_scat_tab_get, _tab_interp.AGN_heat_table_agn_opac_scat_tab_set)
    agn_opac_abs_tab = property(_tab_interp.AGN_heat_table_agn_opac_abs_tab_get, _tab_interp.AGN_heat_table_agn_opac_abs_tab_set)
    agn_arad_tab = property(_tab_interp.AGN_heat_table_agn_arad_tab_get, _tab_interp.AGN_heat_table_agn_arad_tab_set)
    agn_line_co1 = property(_tab_interp.AGN_heat_table_agn_line_co1_get, _tab_interp.AGN_heat_table_agn_line_co1_set)
    agn_line_co2 = property(_tab_interp.AGN_heat_table_agn_line_co2_get, _tab_interp.AGN_heat_table_agn_line_co2_set)
    agn_line_hcn1 = property(_tab_interp.AGN_heat_table_agn_line_hcn1_get, _tab_interp.AGN_heat_table_agn_line_hcn1_set)
    agn_line_hcn2 = property(_tab_interp.AGN_heat_table_agn_line_hcn2_get, _tab_interp.AGN_heat_table_agn_line_hcn2_set)
    agn_line_h2_1 = property(_tab_interp.AGN_heat_table_agn_line_h2_1_get, _tab_interp.AGN_heat_table_agn_line_h2_1_set)
    agn_line_h2_2 = property(_tab_interp.AGN_heat_table_agn_line_h2_2_get, _tab_interp.AGN_heat_table_agn_line_h2_2_set)
    agn_line_h2_3 = property(_tab_interp.AGN_heat_table_agn_line_h2_3_get, _tab_interp.AGN_heat_table_agn_line_h2_3_set)
    agn_line_12m = property(_tab_interp.AGN_heat_table_agn_line_12m_get, _tab_interp.AGN_heat_table_agn_line_12m_set)
    agn_line_8m = property(_tab_interp.AGN_heat_table_agn_line_8m_get, _tab_interp.AGN_heat_table_agn_line_8m_set)
    agn_line_850m = property(_tab_interp.AGN_heat_table_agn_line_850m_get, _tab_interp.AGN_heat_table_agn_line_850m_set)
    agn_column_out_tab = property(_tab_interp.AGN_heat_table_agn_column_out_tab_get, _tab_interp.AGN_heat_table_agn_column_out_tab_set)
    tables = property(_tab_interp.AGN_heat_table_tables_get, _tab_interp.AGN_heat_table_tables_set)
    ntabs = property(_tab_interp.AGN_heat_table_ntabs_get, _tab_interp.AGN_heat_table_ntabs_set)
    lineArrays = property(_tab_interp.AGN_heat_table_lineArrays_get, _tab_interp.AGN_heat_table_lineArrays_set)
    agn_ntemp = property(_tab_interp.AGN_heat_table_agn_ntemp_get, _tab_interp.AGN_heat_table_agn_ntemp_set)
    agn_ndense = property(_tab_interp.AGN_heat_table_agn_ndense_get, _tab_interp.AGN_heat_table_agn_ndense_set)
    agn_nintensity = property(_tab_interp.AGN_heat_table_agn_nintensity_get, _tab_interp.AGN_heat_table_agn_nintensity_set)
    agn_ncolumn_in = property(_tab_interp.AGN_heat_table_agn_ncolumn_in_get, _tab_interp.AGN_heat_table_agn_ncolumn_in_set)
    agn_temp_vals = property(_tab_interp.AGN_heat_table_agn_temp_vals_get, _tab_interp.AGN_heat_table_agn_temp_vals_set)
    agn_dense_vals = property(_tab_interp.AGN_heat_table_agn_dense_vals_get, _tab_interp.AGN_heat_table_agn_dense_vals_set)
    agn_intensity_vals = property(_tab_interp.AGN_heat_table_agn_intensity_vals_get, _tab_interp.AGN_heat_table_agn_intensity_vals_set)
    agn_column_in_vals = property(_tab_interp.AGN_heat_table_agn_column_in_vals_get, _tab_interp.AGN_heat_table_agn_column_in_vals_set)

    def setupTable(self, labelFile, tableFile, convertLines):
        return _tab_interp.AGN_heat_table_setupTable(self, labelFile, tableFile, convertLines)

    def agn_tab_index(self, id, it, ii, _is):
        return _tab_interp.AGN_heat_table_agn_tab_index(self, id, it, ii, _is)

    def __init__(self):
        _tab_interp.AGN_heat_table_swiginit(self, _tab_interp.new_AGN_heat_table())
    __swig_destroy__ = _tab_interp.delete_AGN_heat_table

# Register AGN_heat_table in _tab_interp:
_tab_interp.AGN_heat_table_swigregister(AGN_heat_table)

class CoolHeatTab(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def interpTab(self, density, temperature, intensity, column_in):
        return _tab_interp.CoolHeatTab_interpTab(self, density, temperature, intensity, column_in)

    def interpTabArray(self, n0, n1, n2, n3):
        return _tab_interp.CoolHeatTab_interpTabArray(self, n0, n1, n2, n3)

    def __init__(self, *args):
        _tab_interp.CoolHeatTab_swiginit(self, _tab_interp.new_CoolHeatTab(*args))
    __swig_destroy__ = _tab_interp.delete_CoolHeatTab

# Register CoolHeatTab in _tab_interp:
_tab_interp.CoolHeatTab_swigregister(CoolHeatTab)



