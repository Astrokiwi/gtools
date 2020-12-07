import os
import json

this_dir, this_filename = os.path.split(__file__)

def setup_dirs() :
    config_path = os.path.join(this_dir, '../prams/directories.json')
    config_dict = dict()
    print("Directories have not been configured for this install, please enter them now.")
    print("You can run this subroutine and re-enter these directories at any time with `python -m gtools.config`")
    print("1. Please enter the directory for movie files to be output")
    config_dict["movie"] = input()
    print("2. Please enter the directory for Cloudy-based dust tables")
    config_dict["table"] = input()
    print("3. Please enter the directory or directories where GIZMO/Gadget runs are stored. Separate multiple directories with a comma `,`")
    x = input()
    if ',' in x:
        config_dict['gizmo'] = ','.split(x)
    else:
        config_dict['gizmo'] = x

    with open(config_path, 'w') as f :
        json.dump(config_dict, f)

    return


if __name__ == '__main__' :
    setup_dirs()
