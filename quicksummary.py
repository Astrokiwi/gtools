import sys
import numpy as np
import h5py
import gizmo_tools
from multiprocessing import Pool
import itertools as it

fields_to_investigate = ("KineticEnergy","InternalEnergy","Density","RadiativeAcceleration","SmoothingLength","TimeStep","cs_p","vcirc_p")

def get_derived_data(f,field):
    if field=="KineticEnergy":
        return np.sum(np.array(f["/PartType0/Velocities"])**2,axis=1)*np.array(f["/PartType0/Masses"])/2*1.9891e53 # convert to erg
#         return np.sum(np.array(f["/PartType0/Velocities"])**2,axis=1)*1.9891e53*.5*1.e-14 # convert to erg
#         return np.sum(np.array(f["/PartType0/Velocities"])**2,axis=1) # km**2/s**2
    elif field=="cs_p":
        return np.sqrt(np.array(f["/PartType0/InternalEnergy"])*1.e10)
    elif field=="vcirc_p":
        return np.sqrt(np.sum(np.array(f["/PartType0/Velocities"])**2,axis=1))*1.e5
    else:
        raise ValueError("Derived value {} not known".format(key))

def summarise_cool_bits(run_id,output_dir,snap_id):
    gizmoDir = gizmo_tools.getGizmoDir(run_id)
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    with h5py.File(fullDir+"/snapshot_"+snap_id+".hdf5","r") as f:

        header = f["/Header"]
        time = header.attrs.get("Time")
        time*= 0.9778e9 # to yr
        time/=1.e6 # to Myr
        
        output = [time,f["/PartType0/Density"].size]
        for field in fields_to_investigate:
            key = "/PartType0/"+field
            if key in f:
                x = np.array(f["/PartType0/"+field])
            else:
                x = get_derived_data(f,field)
            for operation in np.min,np.mean,np.median,np.max:
                output.append(operation(x))
        return output

def summarise_and_dump(run_id,output_dir,dumpsOrdered=True):
    print("Dumping full evolution",run_id,output_dir)

    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id,output_dir,dumpsOrdered=dumpsOrdered)
#     snapf = 100
    snap_ids = ["%03d"%x for x in range(snapf)]

    with Pool(processes=80) as pool:
        output_data = pool.starmap(summarise_cool_bits,zip(it.repeat(run_id),it.repeat(output_dir),snap_ids))
    output_data = np.array(output_data)
    header = "1time 2N"
    icol = 2
    summary_types = ("min","mean","median","max")
    for field in fields_to_investigate:
        for summary_type in summary_types:
            icol+=1
            header+=str(icol)+field+"_"+summary_type+" "
    np.savetxt("data/quicksummary"+run_id+output_dir+".dat",output_data,header=header)


if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]

    summarise_and_dump(run_id,output_dir)