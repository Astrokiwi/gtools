from ..tools import gizmo_tools

run_names = [
"longrun_medflow_settled_defaultaniso_polar",
"longrun_medflow_vesc_defaultaniso",
"longrun_medflow_vesc_defaultaniso_polar",
"longrun_weakflow_rapid_defaultaniso",
"longrun_weakflow_rapid_defaultaniso_polar",
"longrun_weakflow_settled_defaultaniso",
"longrun_weakflow_settled_defaultaniso_polar",
"longrun_weakflow_vesc_defaultaniso",
"longrun_weakflow_vesc_defaultaniso_polar",
"newflow_settled_thin_up",
"newflow_vesc_thin_45",
"newflow_vesc_thin_side",
"newflow_vesc_thin_up"]


def setup(dir=""):
    global run_parameters,ordered_keys
    
    run_parameters = gizmo_tools.load_run_parameters("3032",dir=dir)

    run_parameters = {x:y for x,y in run_parameters.items() if x in run_names}

    ordered_keys = gizmo_tools.run_parameters_names(run_parameters)

def run_names_bigstring():
    print("for i in {} ; do echo $i; done".format(" ".join(run_names)))