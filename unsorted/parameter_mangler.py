import sys
import itertools
import copy
import os

"""parameter_mangler.py

This code takes a sample parameter file, a "mangler" input file, and creates permutations of parameter files based on that. It then outputs the command line commands to execute the parameters
"""

class manglon:
    def __init__(self,first_line,last_line=None):
        self.mangle_texts = dict()
        self.nmangle = 0
        self.first_line = first_line
        if last_line is None:
            self.last_line = first_line
        else:
            self.last_line = last_line
        self.nlines = self.last_line-self.first_line+1
        if self.nlines<=0:
            raise ValueError("first_line>=last_line in manglon")
    
    def add_text(self,label,text):
        self.nmangle+=1
        self.mangle_texts[label]=text

def load_manglia(mangle_file):
    manglia = []
    
    # build mangle structure
    inlines = open(mangle_file).readlines()
    inlines = [line.strip() for line in inlines]
    iline = 0
    curent_manglon = None
    while iline<len(inlines):
        line = inlines[iline]
#         print(iline,line)
        if "-" in line:
            a,b = line.split("-")
            current_manglon = manglon(int(a)-1,int(b)-1)
        else:
            current_manglon = manglon(int(line)-1)
        if current_manglon.nlines==1:
            iline+=1
            while len(inlines[iline])>0:
                line = inlines[iline]
                label,text = line.split(":",1)
                current_manglon.add_text(label,text)
                iline+=1
        else:
            iline+=1
            while iline<len(inlines) and len(inlines[iline])>0:
                line = inlines[iline]
                label,text0 = line.split(":",1)
                text_rest = inlines[iline+1:iline+current_manglon.nlines]
                text = [text0]+text_rest
                current_manglon.add_text(label,text)
                iline+=current_manglon.nlines
        manglia+=[current_manglon]
        iline+=1
    return manglia

def mangle_combinations(manglia):
    key_bases = [x.mangle_texts.keys() for x in manglia]
    key_combinations = list(itertools.product(*key_bases))
    return key_combinations

def filebase_from_combination(output_prefix,key_combination):
    return "_".join([output_prefix]+list(key_combination))

def mangle_pram_files(default_pramfile,manglia,nprocs=16,nmaxprocs=172):
    key_combinations = mangle_combinations(manglia)

    n_run_scripts = nmaxprocs//nprocs
    
    default_pramfile = open(default_pramfile).readlines()
    
    iscript=0
    script_lines = [""]*n_run_scripts
    for key_combination in key_combinations:
        key_combination = list(key_combination)
        current_pramfile = copy.deepcopy(default_pramfile)
        output_file = filebase_from_combination(output_prefix,key_combination)#+output_suffix
        for inkey,current_manglon in zip(key_combination,manglia):
            text_choice = current_manglon.mangle_texts[inkey]
            if current_manglon.nlines==1:
                current_pramfile[current_manglon.first_line] = text_choice+"\n"
            else:
                current_pramfile[current_manglon.first_line:current_manglon.last_line+1] = [x+"\n" for x in text_choice]
        current_pramfile[14] = "OutputDir   "+output_file+"\n"
        with open(output_file+".param",'w') as f:
            f.write("".join(current_pramfile))
        script_lines[iscript] += "nice -4 mpirun -np {} ./GIZMO {}.param >& {}.out\n".format(nprocs,output_file,output_file)
        iscript+=1
        iscript=iscript%n_run_scripts
    return n_run_scripts
    
if __name__=='__main__':
    default_pramfile = sys.argv[1]
    mangle_file = sys.argv[2]
    
    output_prefix = "testflows"
    output_suffix = ".param"
        
    manglia = load_manglia(mangle_file)
    
    script_lines = mangle_pram_files(default_pramfile,manglia)

    for iscript,script_runs in enumerate(script_lines):
        with open("runscript_{}.sh".format(iscript),'w') as f:
            f.write(script_runs)



