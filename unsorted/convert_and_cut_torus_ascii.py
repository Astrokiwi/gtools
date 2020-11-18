import gizmo_tools
import sys
import numpy as np

if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    snap_str = sys.argv[3]

    t,particles = gizmo_tools.load_gizmo_pandas(
                                run_id,output_dir,snap_str,
                                ["Masses","Coordinates","Velocities","SmoothingLength","InternalEnergy"],
                                internal_units = True)
    
    rad3d = np.sqrt(particles["Coordinates_x"]**2\
              +particles["Coordinates_y"]**2\
              +particles["Coordinates_z"]**2)

    particles=particles[rad3d>particles["SmoothingLength"]].reset_index(drop=True) 
    
#     particles = particles[:100]
    
    gizmo_tools.dump_ascii(
               f"data/{run_id}_{output_dir}_{snap_str}.txt",
               particles)  
