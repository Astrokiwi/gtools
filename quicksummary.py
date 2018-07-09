import sys
import numpy as np
import h5py
import gizmo_tools
from multiprocessing import Pool
import itertools as it

fields_to_investigate = ("InternalEnergy","Density","RadiativeAcceleration","SmoothingLength","TimeStep")

def summarise_cool_bits(run_id,output_dir,snap_id):
    gizmoDir = gizmo_tools.getGizmoDir()
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    with h5py.File(fullDir+"/snapshot_"+snap_id+".hdf5","r") as f:

        header = f["/Header"]
        time = header.attrs.get("Time")
        time*= 0.9778e9 # to yr
        time/=1.e6 # to Myr
        
        output = [time]
        for field in fields_to_investigate:
            x = np.array(f["/PartType0/"+field])
            for operation in np.min,np.mean,np.median,np.max:
                output.append(operation(x))
        return output


if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]

    print("Dumping full evolution",run_id,output_dir)

    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir)
    snap_ids = ["%03d"%x for x in range(snapf)]

    with Pool(processes=80) as pool:
        output_data = pool.starmap(summarise_cool_bits,zip(it.repeat(run_id),it.repeat(output_dir),snap_ids))
    output_data = np.array(output_data)
    header = "1time "
    icol = 1
    summary_types = ("min","mean","median","max")
    for field in fields_to_investigate:
        for summary_type in summary_types:
            icol+=1
            header+=str(icol)+field+"_"+summary_type+" "
#     print(output_data)
    np.savetxt("data/quicksummary"+run_id+output_dir+".dat",output_data,header=header)
