"""

"""

# External library
import subprocess
import numpy as np
import scipy.spatial as spatial
import MDAnalysis as mda
from numba import jit
from math import floor
import networkx as nx
import scipy.spatial as spt


# Internal library
import initialisation
import thermalisation
import simulation

class ProgPoly():
    # Take the path to gmx, and the initialisation/thermalisation/simulation objects
    def __init__(this, gromacs, initialisation, thermalisation, simulation, ncore, offset, iswarning = True, isverbose = True, isgmxverbose = False):
        this.warning = iswarning # Activate/Desactivate warnings
        this.verbose = isverbose # Activate/Desactivate informative message
        this.gmxverbose = isgmxverbose # Make gmx stop talking (loud and noisy according to their own doc https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html)
        this.gmx = gromacs # gmx app path
        this.I = initialisation # Initialisation object
        this.T = thermalisation # Thermalisation object
        this.S = simulation # Simulation object
        this.ncore = ncore # Number of allocated core
        this.offset = offset # Pin offset
        this.__is_stopped = False # Nothing will work when true
        this.__gmx_working = False # Nothing will work unless true
        this.folder="_DEFAULT_" # Unique folder in which all the intermediate file will be stored

    # If warning is True, print what
    def __warn(this, what):
        if this.warning:
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nPROG_POLY WARNING: "+what+"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    # If verbose is True, print what
    def __speak(this, what):
        if this.verbose:
            print("PROG_POLY "+this.folder+" : "+what)
    # Used to stop the program, no further call of gmx will work
    def __stop(this, error):
        print("PROG_POLY: The program had to stop due to a fatal error: " + error)
        this.__is_stopped=True
    # Used to call a command
    def call(this, command, critical = True):
        if this.__gmx_working and not this.__is_stopped:
            this.__speak("Executing command : "+command+" ...")
            try:
                if this.gmxverbose:
                    subprocess.run(this.gmx + " -nobackup -quiet " + command[3:], shell=True, check=True)
                else:
                    subprocess.run(this.gmx + " -nobackup -quiet " + command[3:], shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                this.__speak("... succeed")
                return True
            except subprocess.CalledProcessError:
                this.__speak("... failed")
                if critical:
                    this.__stop("gmx error when using the command : " + command)
                return False

    # Try/Catch to test if gmx can be launch without problems
    def test(this):
        if not this.__gmx_working:
            this.__speak("Trying to reach gmx")
            try:
                subprocess.run(this.gmx + " -nobackup -quiet", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                this.__speak("gmx found")
                this.__gmx_working = True
            except subprocess.CalledProcessError:
                this.__stop("gmx not found")
                return False
        ps = this.I.default_parameters
        Np = ps["number_p"][0]
        Nm = ps["number_m"][0]
        qp = ps["charge_p"][0]
        qm = ps["charge_m"][0]
        this.__speak("Testing the system charge")
        Q = Np*qp+Nm*qm
        if Q == 0:
            this.__speak("The system has no charge")
        else:
            this.__stop("The system has a charged of " + str(Q) + "e")
            return False
        this.__speak("Checking if the name is correct")
        if this.folder=="_DEFAULT_":
            this.__stop("No name specified, please use the .name(dictionnary) function to specified one. For example .name(q=2, E=1) will make the simulation in a folder 'q_2E_1' and will save the output with similar name")
            return False
        else:
            this.__speak("Correct name " + this.folder)
        return True
    # Function to change some parameter inside file (I will use the form {old} -> new for each old in olds and new in news)
    def change(this, template, file, olds, news):
        lines = ""
        with open(template, "r") as f:
            # Read the lines
            lines = f.read()
        # Finding and replacing
        for i in range(len(olds)):
            lines = lines.replace("{"+str(olds[i])+"}", str(news[i]))
        # Rewrite the file
        with open(file, "w") as f:
            f.write(lines)
    # Main function of the class, launch all the simulation
    def launch(this):
        subprocess.run("mkdir "+this.folder, shell=True)
        this.__speak("Starting the initialisation step")
        # Testing if gmx is accessible
        if not this.test():
            return False
        # Initialisation step
        if not this.I.__launch(this):
            return False
        this.__speak("The system is ready for thermalisation!")
        if not this.T.__launch(this):
            return False
        this.__speak("The system is ready for the simulation!")
        if not this.S.__launch(this):
            return False
        this.__speak("The simulation is finished!")
    def launch_init(this):
        subprocess.run("mkdir "+this.folder, shell=True)
        this.__speak("Starting the initialisation step")
        # Testing if gmx is accessible
        if not this.test():
            return False
        # Initialisation step
        if not this.I.__launch(this):
            return False
    # Remove all the created file to obtain a fresh folder
    def clean(this):
        subprocess.run("rm -r "+this.folder, shell=True)
        subprocess.run("rm step*", shell=True)
        subprocess.run("rm mdout.mdp", shell=True)
    def name(this, **parameters):
        this.folder=""
        for parameter in parameters:
            this.folder+=parameter+"_"+str(parameters[parameter])
        this.folder+="/"
    def compute_intensity(this, initial_step, folder_out="res/int/"):
        freq = this.S.default_parameters["frequency_output_v"][0]
        i_s = int(initial_step/freq)
        file_in = this.folder + "run"
        file_out = folder_out + this.folder[:-1] + ".txt"
        with open(file_out, "w+") as f:
            f.write("# I (A/ps)\n")
        u = mda.Universe(file_in+".gro", file_in+".trr")
        N_t = len(u.trajectory)
        name_p = "Na"
        name_n = "Cl"
        ion_p = u.select_atoms("name "+name_p)
        ion_n = u.select_atoms("name "+name_n)
        this.__speak("Saving trajectory from step " + str(initial_step) + " to " + str(N_t*freq))
        for t in range(i_s, N_t):
            u.trajectory[t]
            cp = u.select_atoms("name "+name_p)
            vp = cp.velocities[:,0]
            cm = u.select_atoms("name "+name_n)
            vm = cm.velocities[:,0]
            with open(file_out, "a+") as f:
                f.write(str((vp-vm).sum())+"\n")
    def compute_cluster(this, folder_out="res/clu/"):
        #ATTENTION FAUT CHANGER LX/LY SELON LES BESOINS
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

        file_in = this.folder + "run"
        file_out = folder_out + this.folder[:-1] + ".txt"
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
            
        

"""
init = initialisation.Initialisation(True, False)

thermalise = thermalisation.Thermalisation(True, True)

run = simulation.Simulation(True, True)
run.fix(number_step=10000)

progpoly = ProgPoly("~/Software/gromacs/build/bin/gmx", init, thermalise, run, 8, 0, True, True, False)
progpoly.name(q=1, E=0)
progpoly.launch()
progpoly.compute_intensity(3000)
progpoly.clean()
"""