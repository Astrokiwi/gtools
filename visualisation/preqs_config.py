lines = ["co1","co2","hcn1","hcn2","h2_1","h2_2","h2_3","12mic","8mic","850mic"]

preqs_list = [
     [["vmag","vel_2d","vel_r","vel_x","vel_y","vel_z","vel_a"]+["v"+line for line in lines]+["vels"+line for line in lines],"vels"]
    ,[["view" ,"dusttau" ,"emit"] , ["tdust","opac","dg"]]
    ,["facetemp" , ["tdust","opac"]]
    ,["vlos"  , ["opac","dg"]]
    ,[["vlos" ,"vmag" ,"vthin"] , "vels"]
    ,["emit"  , ["tdust","dg"]]
    ,["tdust"  , "table"]
    ,["dust"  , "dg"]
    ,[["IRdust","IRdustm"], "dg"]
    ,["dg"  , "table"]
    ,["table"  ,  ["temp","nH","tau"]]
    ,[["vel0" ,"age"],"id"]
    ,["arad","arads"]
    ,[lines,"table"]
]

for line in lines:
    preqs_list+=[["dv"+line,"v"+line]]
    preqs_list+=[["v"+line,["view"+line,"opac","dg"]]]
    preqs_list+=[["view"+line,[line+"m","opac","dg"]]]
    preqs_list+=[["vels"+line,[line+"m","vels"]]]
    preqs_list+=[[line+"m",line]]

preqs_dict = dict()

def add_key(key,preqs_dict):
    if key in preqs_dict: # append if already there
        old_preqs = preqs_dict[key]
        if isinstance(old_preqs,list):
            if isinstance(preqs,list):
                preqs_dict[key]+=preqs
            else:
                preqs_dict[key]+=[preqs]
        else:
            if isinstance(preqs,list):
                preqs_dict[key]=[old_preqs]+preqs
            else:
                preqs_dict[key]=[old_preqs,preqs]
    else: # otherwise add it
        preqs_dict[key]=preqs

for row in preqs_list:
    keys,preqs = row
    if isinstance(keys,list):
        for key in keys:
            add_key(key,preqs_dict)
    else:
        add_key(keys,preqs_dict)


def process_preqs(need_to_load):
    """process_preqs(need_to_load)

    Adds prerequisites to input list or set of strings of variables, returns set including all prerequisites"""

    n_to_load = 0
    while n_to_load!=len(need_to_load): # i.e. unchanged
        n_to_load=len(need_to_load)
        new_need_to_load = need_to_load.copy()
        for to_load in need_to_load:
            if to_load in preqs_dict:
                all_to_add = preqs_dict[to_load]
                if isinstance(all_to_add,list):
                    new_need_to_load.update(all_to_add)
                else:
                    new_need_to_load.add(all_to_add)
        need_to_load = new_need_to_load

    return need_to_load


