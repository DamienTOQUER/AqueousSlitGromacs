"""


"""
# External library
import subprocess

# Old internal library, they come from the old program (they never call gmx, they only create new file)
import gen_slit
import hcp

class Initialisation():
    # Units
    U = {
        "length" : "nm",
        "charge" : "e",
        "none" : ""
    }
    # Constructor, set if the warning should be print (default True) and if the verbose should be on (default False)
    def __init__(this, iswarning = True, isverbose = False):
        this.warning = iswarning
        this.verbose = isverbose
        # Default parameters for the initialisation step, each parameter is an array where the first data is its value and the second its units.
        U = this.U
        this.default_parameters = {
        # Length parameters in nm
        "length_x" : [20,U["length"]],
        "length_y" : [20,U["length"]],
        "length_z" : [20,U["length"]],
        "height" : [0.7,U["length"]],
        "spacing_cc" : [0.142,U["length"]],
        "spacing_hcp" : [0.3,U["length"]],
        # Charge in e
        "charge_p" : [1,U["charge"]],
        "charge_m" : [-1,U["charge"]],
        # Number
        "number_p" : [100,U["none"]],
        "number_m" : [100,U["none"]],
        "number_w" : [4500,U["none"]],
        "SEED" : [0,U["none"]]
        }
        this.ion_structure = {
        # Boolean
        "excluded" : [False, U["none"]],
        "freezed" : [False, U["none"]],
        # Length parameters in nm
        "positions_p" : [[], U["length"]],
        "positions_m" : [[], U["length"]]
        }

    # Private function which change one parameters "name" to a new "value"
    def __replace(this, name, value):
        if name in this.default_parameters:
            a = this.default_parameters[name][0]
            this.default_parameters[name][0] = value
            this.__speak("parameter " + name + " changed from " + str(a) + str(this.default_parameters[name][1]) + " to " + str(this.default_parameters[name][0]) + str(this.default_parameters[name][1]))
        else:
            this.__warn("unknown parameter: " + name)

    # Public function to change multiple parameters (dictionnary)
    def fix(this, **parameters):
        for name in parameters:
            this.__replace(name, parameters[name])

    # If warning is True, print what
    def __warn(this, what):
        if(this.warning):
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nPROG_POLY INITIALISATION WARNING: "+what+"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

    # If verbose is True, print what
    def __speak(this, what):
        if(this.verbose):
            print("PROG_POLY INITIALISATION: "+what)

    # Private function which print the parameter "name" with its unit
    def __show(this, name):
        print(name + ": " + str(this.default_parameters[name][0]) + str(this.default_parameters[name][1]))

    # Public function to show multiple parameters or group of parameter
    def show(this, *names):
        ps = this.default_parameters
        for name in names:
            # For a parameter directly
            if name in ps:
                this.__speak("Showing parameters corresponding to " + name)
                this.__show(name)
            # Show all information
            elif name == "all":
                this.__speak("Showing parameters corresponding to " + name)
                for parameter in ps:
                    this.__show(parameter)
            # Show information according to a certain unit
            elif name in this.U:
                this.__speak("Showing parameters corresponding to " + name)
                for parameter in ps:
                    if ps[parameter][1] == this.U[name]:
                        this.__show(parameter)
            # Show information about the simulation box
            elif name == "box":
                this.__speak("Showing parameters corresponding to the simulation box")
                lx = ps["length_x"]
                ly = ps["length_y"]
                lz = ps["length_z"]
                print("The box dimension is " + str(lx[0]) + str(lx[1]) + "*" + str(ly[0]) + str(ly[1]) + "*" + str(lz[0]) + str(lz[1]) + " for a total volume of " + str(lx[0]*ly[0]*lz[0]) + str(lx[1])+"^3")
            # Show information about the slit
            elif name == "slit":
                this.__speak("Showing parameters corresponding to the slit")
                lx = ps["length_x"]
                ly = ps["length_y"]
                lz = ps["height"]
                dcc = ps["spacing_cc"]
                nw = ps["number_w"]
                print("The slit dimension is " + str(lx[0]) + str(lx[1]) + "*" + str(ly[0]) + str(ly[1]) + "*" + str(lz[0]) + str(lz[1]) + " for a total volume of " + str(lx[0]*ly[0]*lz[0]) + str(lx[1])+"^3. It is composed of 2 graphene sheets, with a distance between carbon atoms of " + str(dcc[0]) + str(dcc[1]) + ". There is " + str(nw[0]) + str(nw[1]) + " water molecules, which correpond to a density of " + str(nw[0]/lx[0]/ly[0]/lz[0]) + str(lx[1]) +"^-3")
            elif name == "numerical":
                this.__speak("Showing parameters corresponding to numerical aspects of the simulation")
                lp = ps["spacing_hcp"]
                print("The seed is " + str(ps["SEED"][0]) + ". The hcp lattice (used to place water and ions) has a lattice parameter of " + str(lp[0]) + str(lp[1]))
            elif name == "ions":
                this.__speak("Showing parameters corresponding to the ions")
                cp = ps["charge_p"]
                cm = ps["charge_m"]
                np = ps["number_p"]
                nm = ps["number_m"]
                print("The cation has a charge of " + str(cp[0]) + str(cp[1]) + ", there is " + str(np[0]) + " of them in the simulation. The anion has a charge of " + str(cm[0]) + str(cm[1]) + ", there is " + str(nm[0]) + " of them in the simulation. The total charge is " + str(np[0]*cp[0]+nm[0]*cm[0]) + str(cp[1]))
    # Public function to show the units system used in the simulation
    def units(this):
        print("The unit system is the following:")
        for u in this.U:
            if u!="none":
                print("\t" + u + " is expressed in " + this.U[u])
    # Public function which add ions, compute the center of mass and the max radius
    def add_ion_structure(this, pos_p, pos_m, isfreezed=True):
        this.ion_structure["excluded"][0] = True
        this.ion_structure["positions_p"][0] = pos_p
        this.ion_structure["positions_m"][0] = pos_m
        this.ion_structure["freezed"][0] = isfreezed

    # Main function, used only from the progpoly class "pp"
    def _ProgPoly__launch(this, pp):
        ps = this.default_parameters
        # Getting the needed parameters
        d = ps["spacing_cc"][0]
        h = ps["height"][0]
        Lx = ps["length_x"][0]
        Ly = ps["length_y"][0]
        Lz = ps["length_z"][0]
        d_w = ps["spacing_hcp"][0]
        Nw = ps["number_w"][0]
        Np = ps["number_p"][0]
        Nm = ps["number_m"][0]
        SEED = ps["SEED"][0]

        if this.ion_structure["excluded"][0]:
            Np += this.ion_structure["positions_p"][0].shape[0]
            Nm += this.ion_structure["positions_m"][0].shape[0]

        ## Generating the graphene slit
        # Create the motif.gro file
        this.__speak("Generation of the initial configuration")
        dx, dy, dz = gen_slit.gen_motif(d, h, pp.folder+"motif.gro")
        this.__speak("Calculating the number of motif along each direction to obtain a good length")
        Nx = int(Lx/dx) # Number of motif along x
        Ny = int(Ly/dy) # Number of motif along y
        # We adjust the length
        """ DO NOT DO THAT, DUE TO NUMERICAL IMPRECISION THE SYSTEM WILL BECOME SMALLER AND SMALLER IF WE INITIALISE IT AGAIN
        ps["length_x"][0] = Nx*dx
        ps["length_y"][0] = Ny*dy
        Lx = ps["length_x"][0]
        Ly = ps["length_y"][0]
        """
        lx = Nx*dx
        ly = Ny*dy
        # We then generate the slit using genconf
        pp.call('gmx genconf -f '+pp.folder+'motif.gro -o '+pp.folder+'slit.gro -nbox ' + str(Nx) + ' ' + str(Ny) + ' 1')
        # Change the box size in the file (to have a correct length along z)
        with open(pp.folder+"slit.gro", "r+") as f:
            # Read the lines
            lines = f.readlines()
            # Begining of the file
            f.seek(0)
            # Empty the file
            f.truncate()
            # Remove the box volume at the end
            last = lines[-1]
            # Change the box volume
            lastab = last.split(" ")
            lastab[-1] = str(Lz)+"\n"
            nothing = " " #It's nothing
            lines[-1] = nothing.join(lastab)
            # Rewrite the file
            f.writelines(lines)
        this.__speak("The graphene slit has been created")

        ## Solvation and ionisation of the system
        # Calculation of the center of mass and max radius of the ion structure
        pos_p = None
        pos_m = None
        pos_com = None
        radius = None
        if this.ion_structure["excluded"][0]:
            pos_p = this.ion_structure["positions_p"][0]
            pos_m = this.ion_structure["positions_m"][0]
            pos_com = (pos_p.sum(axis=0)+pos_m.sum(axis=0))/(pos_p.shape[0]+pos_m.shape[0])
            radius = max(((pos_p-pos_com)**2).sum(axis=1).max(), ((pos_m-pos_com)**2).sum(axis=1).max())**0.5+d_w
        # Generation of the hcp lattice
        this.__speak("Creation of the hcp lattice")
        tab_w = None
        tab_p = None
        tab_m = None
        if not this.ion_structure["excluded"][0]:
            tab_w, tab_p, tab_m = hcp.gen_hcp_cz(d_w, d, d, d, lx-d, ly-d, h+d, Nw, Np, Nm, SEED)
        else:
            tab_w, tab_p, tab_m = hcp.gen_hcp_cz_excludecircle(d_w, d, d, d, lx-d, ly-d, h+d, Nw, Np-this.ion_structure["positions_p"][0].shape[0], Nm-this.ion_structure["positions_m"][0].shape[0], pos_com[0], pos_com[1], radius, SEED)
            tab_p = pos_p
            tab_m = pos_m
            this.__warn("For the moment, if ion are added using add_ion_structure, only these ions can exist inside the system.")
        # Solvation
        if Nw>0:
            gen_slit.add_water(tab_w, SEED, pp.folder+"slit.gro")
            this.__speak("The system has been hydrated")
        else:
            this.__warn("No water found")
        # Adding the ions
        gen_slit.add_ions(tab_p, pp.folder+"slit.gro", "IP", "Na")
        gen_slit.add_ions(tab_m, pp.folder+"slit.gro", "IN", "Cl")
        # Changing the ff.itp file
        pp.change("template/ff.itp", pp.folder+"ff.itp", ["charge_p", "charge_m"], [ps["charge_p"][0], ps["charge_m"][0]])
        this.__speak("The system has been ionised")

        ## Generate the topology file
        with open(pp.folder+"topol.top", "w") as f:
            f.write('#include "ff.itp"\n')
            f.write('#include "../itp/ions.itp"\n')
            f.write('#include "../itp/carbons.itp"\n')
            f.write('#include "../itp/spce.itp"\n')
            f.write('[ system ]\n')
            f.write('Water molecule and NaCl between 2 graphene sheets\n')
            f.write('[ molecules ]\n')
            f.write('CB '+str(8*Nx*Ny)+'\n')
            f.write('SOL '+str(Nw)+'\n')
            f.write('IP '+str(Np)+'\n')
            f.write('IN '+str(Nm)+'\n')
        this.__speak("The topology file has been generated")

        ## Generate the index .ndx file
        Ntot = 8*Nx*Ny+Nw*3+Np+Nm
        nothing = "\n"
        System = nothing.join([str(i+1) for i in range(Ntot)])
        CB = nothing.join([str(i+1) for i in range(8*Nx*Ny)])
        SOL = nothing.join([str(i+1) for i in range(8*Nx*Ny, 8*Nx*Ny+3*Nw)])
        IONS = nothing.join([str(i+1) for i in range(8*Nx*Ny+3*Nw, Ntot)])
        IP = nothing.join([str(i+1) for i in range(8*Nx*Ny+3*Nw, 8*Nx*Ny+3*Nw+Np)])
        IM = nothing.join([str(i+1) for i in range(8*Nx*Ny+3*Nw+Np, Ntot)])
        Notfreeze = None
        Freeze = None
        if not this.ion_structure["freezed"][0]:
            Notfreeze = nothing.join([str(i+1) for i in range(8*Nx*Ny, Ntot)])
            Freeze = nothing.join([str(i+1) for i in range(8*Nx*Ny)])
        else:
            Notfreeze = nothing.join([str(i+1) for i in range(8*Nx*Ny, 8*Nx*Ny+3*Nw)])
            Freeze = nothing.join([str(i+1) for i in range(8*Nx*Ny)]+[str(i+1) for i in range(8*Nx*Ny+3*Nw, Ntot)])
            this.__warn("When freezed ion are added, there are put inside the Notfreeze group, so if you want them to become dynamic again later, you have to modify the index.ndx groups.")
        res = nothing.join(["[ System ]", System, "[ CB ]", CB, "[ SOL ]", SOL, "[ IONS ]", IONS, "[ Notfreeze ]", Notfreeze, "[ IP ]", IP, "[ IM ]", IM, "[ Freeze ]", Freeze])
        res+=nothing
        with open(pp.folder+"index.ndx", "w") as f:
            f.writelines(res)
        this.__speak("The index file has been generated")

        ## Minimize the structure
        max_try = 10
        this.__speak("Minimalization step")
        other_key = ["freeze_group", "freeze_dim"]
        other_value = None
        other_value = ["CB", "Y Y Y"]
        if not this.ion_structure["freezed"][0]:
            other_value = ["CB", "Y Y Y"]
        else:
            other_value = ["CB IONS", "Y Y Y Y Y Y"]
        pp.change("template/em.mdp", pp.folder+"em.mdp", ["SEED"]+other_key, [ps["SEED"][0]]+other_value)
        pp.call("gmx grompp -f "+pp.folder+"em.mdp -c "+pp.folder+"slit.gro -r "+pp.folder+"slit.gro -p "+pp.folder+"topol.top -n "+pp.folder+"index.ndx -o "+pp.folder+"em.tpr")
        if not pp.call("gmx mdrun -v -deffnm "+pp.folder+"em -nt " + str(pp.ncore) + " -pin on -pinoffset " + str(pp.offset) + " -pinstride 2", False):
            this.__warn("The energy minimisation returned a segment fault (assumption, maybe it's something else). The program will try again up to "+str(max_try)+" times with other seeds")
            i = 0
            while i<max_try-2:
                pp.change("template/em.mdp", pp.folder+"em.mdp", ["SEED"]+other_key, [ps["SEED"][0]+i+1]+other_value)
                pp.call("gmx grompp -f "+pp.folder+"em.mdp -c "+pp.folder+"slit.gro -r "+pp.folder+"slit.gro -p "+pp.folder+"topol.top -n "+pp.folder+"index.ndx -o "+pp.folder+"em.tpr")
                if pp.call("gmx mdrun -v -deffnm "+pp.folder+"em -nt " + str(pp.ncore) + " -pin on -pinoffset " + str(pp.offset) + " -pinstride 2", False):
                    break
                else:
                    i=i+1
                    this.__warn("Try n°" + str(i+1) + " failed")
            if i==max_try-2:
                pp.change("template/em.mdp", pp.folder+"em.mdp", ["SEED"]+other_key, [ps["SEED"][0]+i+1]+other_value)
                pp.call("gmx grompp -f "+pp.folder+"em.mdp -c "+pp.folder+"slit.gro -r "+pp.folder+"slit.gro -p "+pp.folder+"topol.top -n "+pp.folder+"index.ndx -o "+pp.folder+"em.tpr")
                pp.call("gmx mdrun -v -deffnm "+pp.folder+"em -nt " + str(pp.ncore) + " -pin on -pinoffset " + str(pp.offset) + " -pinstride 2")
                return False
            this.__warn("Finally succeed")
        return True

