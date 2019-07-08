import numpy as np
import gizmo_tools
import sys
import os
# import pandas as pd

if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]

    gizmoDir = gizmo_tools.getGizmoDir(run_id)
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    filename = gizmoDir+"/"+run_id+"/"+output_dir+"/cpu.txt"

#     tempFile = "data/cputemp{}{}.dat".format(run_id,output_dir)
#     print("Grepping")
#     os.system("grep Step {} > {}".format(filename,tempFile))
#     print("Reading")
#     time = pd.read_table(filename, sep=',| ', engine='python', header=None)
#     print(time)
#     os.remove(tempFile)

#     sys.exit()
    sim_time = []
    clock_time = []
    iline = 0
    iskip = 100
    with open(filename) as f:
        for line in f:
            iline+=1
            if iline%100000==0: print(iline)
            if "Time" in line:
                a = line.split()
                sim_time.append(float(a[3][:-1]))
            elif "total" in line:
                a = line.split()
                clock_time.append(float(a[1]))
    print(len(sim_time),len(clock_time))
    outp = np.array([sim_time,clock_time]).T
    outp = outp[::iskip,:]
    np.savetxt("data/cputime{}{}.dat".format(run_id,output_dir),np.array([sim_time,clock_time]).T)          