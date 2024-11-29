import numpy as np
import scipy.spatial as spatial

# CREATE THE MOTIF FOR THE BI-GRAPHENE SHEET AND STORE IT IN A .gro FILE
def gen_motif(d, h, output_file):
    l = h
    dx = 3**0.5*d
    dy = 3*d
    with open(output_file, "w") as f:
        # Header
        f.write("Graphene double sheet (slit) motif\n")
        # Atoms
        f.write("8\n")
        # Sheet 1
        f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (1, "CB", "C", 1, 0, 0, d/2, 0, 0, 0))
        f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (1, "CB", "C", 2, 0, d, d/2, 0, 0, 0))
        f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (1, "CB", "C", 3, 3**0.5/2*d, 3/2*d, d/2, 0, 0, 0))
        f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (1, "CB", "C", 4, 3**0.5/2*d, 5/2*d, d/2, 0, 0, 0))
        # Sheet 2
        f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (2, "CB", "C", 5, 0, 0, l+d/2, 0, 0, 0))
        f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (2, "CB", "C", 6, 0, d, l+d/2, 0, 0, 0))
        f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (2, "CB", "C", 7, 3**0.5/2*d, 3/2*d, l+d/2, 0, 0, 0))
        f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (2, "CB", "C", 8, 3**0.5/2*d, 5/2*d, l+d/2, 0, 0, 0))
        # Box size
        f.write("%8.4f%8.4f%8.4f" % (dx, dy, l+d))
    return (dx, dy, l)

# ADD THE WATER MOLECULE (with a random rotation) IN THE SLIT    
def add_water(positions, SEED, file, res_name="SOL", O_name = "OW", H1_name = "HW1", H2_name = "HW2"):
    dOH = 0.1 # nm
    theta = 109.47 # °
    theta = theta/180*np.pi
    Hw1 = np.array([dOH, 0.0, 0.0])
    Hw2 = np.array([dOH*np.cos(theta), dOH*np.sin(theta), 0.0])
    random_rot = spatial.transform.Rotation.random(positions.shape[0], random_state=SEED)
    Hw1_rot = random_rot.as_matrix()@Hw1
    Hw2_rot = random_rot.as_matrix()@Hw2
    with open(file, "r+") as f:
        # Read the lines
        lines = f.readlines()
        # Begining of the file
        f.seek(0)
        # Empty the file
        f.truncate()
        # Remove the box volume at the end
        last = lines[-1]
        # Last residue number
        lm = int(lines[-2][0:5])
        # Last atom number
        la = int(lines[-2][15:20])
        # New number of atoms
        na = int(lines[1])+3*positions.shape[0]
        # Change the number of atoms
        lines[1]=str(na)+"\n"
        # Rewrite the file
        f.writelines(lines[:-1])
        # Add the new water molecules
        for i in range(positions.shape[0]):
            f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (lm+1+i, res_name, O_name, la+1+3*i, positions[i,0], positions[i,1], positions[i,2], 0, 0, 0))
            f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (lm+1+i, res_name, H1_name, la+2+3*i, positions[i,0]+Hw1_rot[i,0], positions[i,1]+Hw1_rot[i,1], positions[i,2]+Hw1_rot[i,2], 0, 0, 0))
            f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (lm+1+i, res_name, H2_name, la+3+3*i, positions[i,0]+Hw2_rot[i,0], positions[i,1]+Hw2_rot[i,1], positions[i,2]+Hw2_rot[i,2], 0, 0, 0))
        # Add the volume at the end
        f.write(last)
    
def add_ions(positions, file, res_name, ion_name):
    with open(file, "r+") as f:
        # Read the lines
        lines = f.readlines()
        # Begining of the file
        f.seek(0)
        # Empty the file
        f.truncate()
        # Remove the box volume at the end
        last = lines[-1]
        # Last residue number
        lm = int(lines[-2][0:5])
        # Last atom number
        la = int(lines[-2][15:20])
        # New number of atoms
        na = int(lines[1])+positions.shape[0]
        # Change the number of atoms
        lines[1]=str(na)+"\n"
        # Rewrite the file
        f.writelines(lines[:-1])
        # Add the new ions
        for i in range(positions.shape[0]):
            f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (lm+1+i, res_name, ion_name, la+1+i, positions[i,0], positions[i,1], positions[i,2], 0, 0, 0))
        # Add the volume at the end
        f.write(last)