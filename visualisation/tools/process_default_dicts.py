import json

def pack_dicts():
    lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3"]
    # count "IRdust" as a line
    lines = ["IRdust"]+lines

    # placeholder definitions
    densslice	=	'densslice'
    weightslice	=	'weightslice'
    zdensslice	=	'zdensslice'
    zweightslice	=	'zweightslice'
    vec2dslice_flat	=	'vec2dslice_flat'
    viewslice	=	'viewslice'
    minslice	=	'minslice'
    zminslice	=	'zminslice'
    vorinoislice	=	'vorinoislice'
    zvorinoislice	=	'zvorinoislice'
    maxslice	=	'maxslice'
    zmaxslice	=	'zmaxslice'
    maxdotslice	=	'maxdotslice'
    mindotslice	=	'mindotslice'
    sdevslice	=	'sdevslice'
    weightviewslice	=	'weightviewslice'
    zvec2dslice	=	'zvec2dslice'


    # "dynamic" properties to deal with later
    quantslice = 'quantslice'
    dslice = 'dslice'
    rhounit = [r"$\log_{10} \Sigma$ (M$_\odot$/pc$^2$)",r"$\log_{10} \rho$ (M$_\odot$/pc$^3$)"]
    drange = [-2.,7.]
    vec2dslice = 'vec2dslice'
    thisminslice = 'thisminslice'
    mslice = 'mslice'

    drange_s = drange


    # Pack the dictionaries!
    # We could do this with a bunch of if/endif statements, but
    # data-driven is better
    plotLabel = dict()
    plotLabel["temp"] = r"$\log_{10} T_g$ (K)"
    plotLabel["col"] = r"$\log_{10} N$ (cm$^{-2}$)"
    plotLabel["nH"] = r"$\log_{10} n_{H}$ (cm$^{-3}$)"
    plotLabel["heat"] = r"$H$"
    plotLabel["dens"] = rhounit
    plotLabel["dusttau"] = r"$\log_{10}\tau_d$"
    plotLabel["depth"] = r"$\tau$"
    plotLabel["vels"] = r"$\log_{10}v$ (km/s)"
    plotLabel["tdust"] = r"$\log_{10} T_d$ (K)"
    plotLabel["dg"] = r"$f_d$"
    plotLabel["dust"] = rhounit
    plotLabel["view"] = r"$\log_{10} F$ (arbitrary units)"
    plotLabel["vlos"] = r"$v_{LOS}$"
    plotLabel["vthin"] = r"$v_{LOS}$"
    plotLabel["emit"] = r"$\log_{10} F$"
    plotLabel["dt"] = r"$\log_{10}$ dt"
    plotLabel["vel_2d"] = r"$v_{2D}$"
    plotLabel["vel_x"] = r"$v_{x}$"
    plotLabel["vel_y"] = r"$v_{y}$"
    plotLabel["vel_z"] = r"$v_{z}$"
    plotLabel["vel_r"] = r"$v_{r}$ (km/s)"
    plotLabel["vmag"] = r"$|v|_{max}$ (km/s)"
    plotLabel["vel_a"] = r"$v_{\theta}$"
    plotLabel["arad"] = r"$a_{rad}$ (log cm/s/s)"
    plotLabel["accel"] = r"$a$ (log cm/s/s)"
    plotLabel["AGNI"] = r"Unextinced $I_{AGN}$ (log erg/s/cm^2)"
    plotLabel["tau"] = r"$\tau$"
    plotLabel["list"] = r"$i$"
    plotLabel["rand"] = r"$q$"
    plotLabel["opac"] = r"$\kappa$"
    plotLabel["smooth"] = r"$h_p$"
    plotLabel["facetemp"] = r"$\log_{10}T_d$ (K)"
    plotLabel["rad0"] = r"$R_0$ (pc)"
    plotLabel["nneigh"] = r"$N_n$"
    plotLabel["pres"] = r"$P$ (internal)"
    plotLabel["age"] = r"$age$ (Myr)"
    for line in lines:
        plotLabel[line]=line
        plotLabel[line+"m"]=line+" line intensity - erg/s/cm**2/ster"
        plotLabel["v"+line]="v"+line
        plotLabel["dv"+line]="dv"+line
        plotLabel["view"+line]="view"+line
        plotLabel["vels"+line]=r"$\log_{10}v_{"+line+r"}$"

    plotRanges = dict()
    plotRanges["temp"] = [.999,6.]*2
    plotRanges["col"] = [18.,26.]*2
    plotRanges["nH"] = [5.,10.]*2
    plotRanges["heat"] = [-1e2,1.e2]*2
    plotRanges["dens"] = drange+drange_s
    plotRanges["dusttau"] = [-3.,3.]*2
    plotRanges["depth"] = [0.,1.e3]*2
    plotRanges["vels"] = [0.,3.]*2
    plotRanges["tdust"] = [0.,3.]*2
    plotRanges["dg"] = [0.0,.000234]*2
    plotRanges["dust"] = drange+drange_s
    plotRanges["view"] = [1.,7.]*2
    plotRanges["vlos"] = [-1.2e2,1.2e2]*2
    plotRanges["vthin"] = [-300.,300.]*2
    plotRanges["emit"] = [-1.,6.]*2
    plotRanges["dt"] = [-3.,3.]*2
    plotRanges["vel_2d"] = [-250.,250.]*2
    plotRanges["vel_r"] = [-300.,300.]*2
    plotRanges["vmag"] = [0.,120.]*2
    plotRanges["vel_x"] = [-200.,200.]*2
    plotRanges["vel_y"] = [-200.,200.]*2
    plotRanges["vel_z"] = [-2000.,2000.]*2
    plotRanges["vel_a"] = [-300.,300.]*2
    plotRanges["arad"] = [-9.,0.]*2
    plotRanges["accel"] = [-3.,0.]*2
    plotRanges["AGNI"] = [-1.,6.]*2
    plotRanges["tau"] = [0.,7.]*2
    plotRanges["list"] = [0.,1.e6]*2
    plotRanges["rand"] = [0.,1.]*2
    plotRanges["opac"] = [-3.,2]*2
    plotRanges["smooth"] = [-4.,2.]*2
    plotRanges["facetemp"] = [0.,5.]*2
    plotRanges["rad0"] = [0.,4.]*2
    plotRanges["nneigh"] = [0,50]*2
    plotRanges["pres"] = [3.3,16]*2
    plotRanges["age"] = [0.,3.]*2
    for line in lines:
        plotRanges[line]=[-8.,2.]*2
        plotRanges[line+"m"]=[-12.,0.]*2
        plotRanges["v"+line]=[-200.,200.]*2
        plotRanges["dv"+line]=[0.,3.,0.,30.]
        plotRanges["view"+line]=[-20.,0.]*2
        plotRanges["vels"+line]=[0.,3.]*2

    plotSliceTypes = dict()
    plotSliceTypes["temp"] = quantslice
    plotSliceTypes["col"] = quantslice
    plotSliceTypes["nH"] = quantslice
    plotSliceTypes["heat"] = quantslice
    plotSliceTypes["dens"] = dslice
    plotSliceTypes["dusttau"] = densslice
    plotSliceTypes["depth"] = zweightslice
    plotSliceTypes["vels"] = vec2dslice
    plotSliceTypes["tdust"] = quantslice
    plotSliceTypes["dg"] = quantslice
    plotSliceTypes["dust"] = dslice
    plotSliceTypes["view"] = viewslice
    plotSliceTypes["vlos"] = viewslice
    plotSliceTypes["emit"] = dslice
    plotSliceTypes["dt"] = thisminslice
    plotSliceTypes["vel_2d"] = quantslice
    plotSliceTypes["vel_r"] = quantslice
    plotSliceTypes["vmag"] = mslice # max
    plotSliceTypes["vel_x"] = quantslice
    plotSliceTypes["vel_y"] = quantslice
    plotSliceTypes["vel_z"] = quantslice
    plotSliceTypes["vthin"] = quantslice
    plotSliceTypes["vel_a"] = quantslice
    plotSliceTypes["arad"] = quantslice
    plotSliceTypes["accel"] = quantslice
    plotSliceTypes["AGNI"] = quantslice
    plotSliceTypes["tau"] = zweightslice
    plotSliceTypes["list"] = quantslice
    plotSliceTypes["opac"] = quantslice
    plotSliceTypes["smooth"] = quantslice
    plotSliceTypes["facetemp"] = viewslice
    plotSliceTypes["rand"] = quantslice
    plotSliceTypes["rad0"] = quantslice
    plotSliceTypes["nneigh"] = quantslice
    plotSliceTypes["pres"] = quantslice
    plotSliceTypes["age"] = quantslice
    for line in lines:
        plotSliceTypes[line]=quantslice
        plotSliceTypes[line+"m"]=dslice
        plotSliceTypes["v"+line]=weightviewslice
        plotSliceTypes["dv"+line]=sdevslice
        plotSliceTypes["view"+line]=viewslice
        plotSliceTypes["vels"+line]=vec2dslice
    
    plotCustomMass = dict()
    plotCustomMass["dust"] = "dust"
    plotCustomMass["dusttau"] = "dusttau"
    for line in lines:
        plotCustomMass[line+"m"]=line+"m"
        plotCustomMass["v"+line]=line+"m"
        plotCustomMass["dv"+line]=line+"m"
        plotCustomMass["vels"+line]=line+"m"
    
    plotData = dict()
    plotData["temp"] = "TK_p"
    plotData["col"] = "coldens"
    plotData["nH"] = "nH_p"
    plotData["heat"] = "heat"
    plotData["depth"] = "depth_p"
    plotData["vels"] = ["vel_x","vel_y","vel2d","vel_z"]
    plotData["tdust"] = "dustTemp"
    plotData["dg"] = "dg"
    plotData["view"] = ["brightness","opac","brightness","opac"]
    plotData["vlos"] = ["vthin","opac","vthin","opac"]
    plotData["emit"] = "emissivity"
    plotData["opac"] = "opac"
    plotData["dt"] = "dt_p"
    plotData["vel_2d"] = "vel2d"
    plotData["vel_r"] = "velr"
    plotData["vmag"] = "vmag"
    plotData["vel_x"] = "vel_x"
    plotData["vel_y"] = "vel_y"
    plotData["vel_z"] = "vel_z"
    plotData["vel_a"] = "vel_a"
    plotData["vthin"] = "vthin"
    plotData["arad"] = "arad"
    plotData["accel"] = "accel"
    plotData["tau"] = "tau"
    plotData["AGNI"] = "AGNI"
    plotData["list"] = "list"
    plotData["rand"] = "rand"
    plotData["smooth"] = "h_p"
    plotData["rad0"] = "rad0"
    plotData["facetemp"] = ["dustTemp","opac","dustTemp","opac"]
    plotData["nneigh"] = "nneigh"
    plotData["pres"] = "pres"
    plotData["age"] = "age"
    for line in lines:
        plotData[line]=line
        plotData[line+"m"]=line+"m"
        plotData["v"+line]=["vthin",line+"brightness",line+"opac","vythin",line+"brightness",line+"opac"]
        plotData["dv"+line]="vthin"
        plotData["view"+line] = [line+"brightness",line+"opac",line+"brightness",line+"opac"]
        plotData["vels"+line]=["vel_x","vel_y","vel2d","vel_z"]
    
    logSliceTypes = [   "temp","col","nH","dens","dust",
                        "view","emit","dt","arad","accel",
                        "AGNI","opac","smooth","tdust","pres","dusttau"]
    logSliceTypes+=lines
    logSliceTypes+=[line+"m" for line in lines]
    logSliceTypes+=["view"+line for line in lines]
    symLogSliceTypes={"heat":1.e-4}

    divergingTypes = ["vel_r","vel_x","vel_y","vel_z","vel_2d","vel_a","vthin"]
    divergingTypes+=["v"+line for line in lines]

    return plotLabel, plotRanges, plotSliceTypes, plotCustomMass, plotData, logSliceTypes, divergingTypes, symLogSliceTypes

plotLabel, plotRanges, plotSliceTypes, plotCustomMass, plotData, logSliceTypes, divergingTypes, symLogSliceTypes = pack_dicts()

field_dict_names = ["label","range","slice","mass","data"]
field_dicts = plotLabel, plotRanges, plotSliceTypes, plotCustomMass, plotData
field_list_names = ["log","diverging","symlog"]
field_lists = logSliceTypes, divergingTypes, symLogSliceTypes

data_keys = field_dicts[0].keys()

big_dict = dict()


for data_key in data_keys:
    big_dict[data_key] = dict()
    for field_dict,field_name in zip(field_dicts,field_dict_names):
        if data_key in field_dict:
            big_dict[data_key][field_name] = field_dict[data_key]
    for field_list,field_name in zip(field_lists,field_list_names):
        if data_key in field_list:
            big_dict[data_key][field_name] = True
        else:
            big_dict[data_key][field_name] = False


with open("plot_defaults.json",'w') as f:
    json.dump(big_dict, f, sort_keys=True, indent=4) 

























