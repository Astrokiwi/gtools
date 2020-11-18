import gizmo_tools
import sys
import numpy as np

if __name__ == '__main__':
    run_id = sys.argv[1]
    output_dir = sys.argv[2]
    snap_str = sys.argv[3]

    t,particles = gizmo_tools.load_gizmo_pandas(
                                run_id,output_dir,snap_str,
                                ["Masses","Coordinates","Velocities","SmoothingLength","InternalEnergy","RadiativeAcceleration"],
                                internal_units = True)
    
    particles["rad3d"] = np.sqrt(particles["Coordinates_x"]**2\
              +particles["Coordinates_y"]**2\
              +particles["Coordinates_z"]**2)

    particles["a3d"] = np.sqrt(particles["RadiativeAcceleration_x"]**2\
              +particles["RadiativeAcceleration_y"]**2\
              +particles["RadiativeAcceleration_z"]**2)

    particles[["rad3d","a3d"]].to_csv(f"data/rad_vs_pressure_{run_id}_{output_dir}_{snap_str}.txt",sep=' ',index=False)