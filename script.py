import subprocess
import numpy as np
import scipy.spatial as spatial
import MDAnalysis as mda
from numba import jit
from math import floor
import networkx as nx
import scipy.spatial as spt

@jit(nopython=True)
def rel_distx(p1, p2):
    Lx = 98.4
    dx = p2[0]-p1[0]
    return dx-floor((dx/Lx+1/2))*Lx

@jit(nopython=True)
def rel_disty(p1, p2):
    Ly = 97.98
    dy = p2[1]-p1[1]
    return dy-floor((dy/Ly+1/2))*Ly

@jit(nopython=True)
def rel_dist(p1, p2):
    return (rel_distx(p1, p2)**2+rel_disty(p1, p2)**2)**0.5

folder = "./q_2.0SEED_1L_20/"
folder_out = "./res/clu/"

file_in = folder + "run"
file_out = folder_out + folder[:-1] + ".txt"
with open(file_out, "w+") as f:
    f.write("Clusters sizes\n")
name_p = "Na"
name_m = "Cl"
u = mda.Universe(file_in+".gro", file_in+".trr")
T = len(u.trajectory)
for t in range(0, T, 10):
    u.trajectory[t]
    ion_p = u.select_atoms("name "+name_p)
    ion_m = u.select_atoms("name "+name_m)
    pos = np.concatenate((ion_p.positions, ion_m.positions), axis=0)
    distance = spt.distance.cdist(pos, pos, rel_dist)
    G = nx.from_numpy_matrix(distance<4.5)
    with open(file_out, "a+") as f:
        for n in range(200):
            f.write(str(len(nx.node_connected_component(G, n)))+" ")
        f.write("\n")