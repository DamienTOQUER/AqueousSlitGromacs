"""


"""
# External library
import subprocess

class Thermalisation():
    # Units
    U = {
        "duration" : "ps",
        "temperature" : "K",
        "field" : "V/nm",
        "length" : "nm",
        "stiffness" : "kJ/mol/nm^2",
        "none" : ""
    }
    # Constructor, set if the warning should be print (default True) and if the verbose should be on (default False)
    def __init__(this, iswarning = True, isverbose = False):
        this.warning = iswarning
        this.verbose = isverbose
        # Default parameters for the initialisation step, each parameter is an array where the first data is its value and the second its units.
        U = this.U
        this.default_parameters = {
        # Duration parameters in ps
        "duration_step" : [0.001,U["duration"]],
        # Temperature in K
        "temperature" : [300,U["temperature"]],
        # Field in V/nm
        "field" : [0, U["field"]],
        # Number
        "number_step" : [10000,U["none"]],
        "frequency_output_x" : [0,U["none"]],
        "frequency_output_v" : [0,U["none"]],
        "frequency_output_f" : [0,U["none"]],
        "frequency_output_e" : [0,U["none"]],
        "frequency_compressed_output" : [0,U["none"]],
        "SEED" : [0,U["none"]]
        }
        this.umbrella = {
        "umbrella" : [False, U["none"]],
        "pair_distance" : [0, U["length"]],
        "stiffness" : [0, U["stiffness"]] 
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
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nPROG_POLY THERMALISATION WARNING: "+what+"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

    # If verbose is True, print what
    def __speak(this, what):
        if(this.verbose):
            print("PROG_POLY THERMALISATION: "+what)

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
    # Public function to show the units system used in the simulation
    def units(this):
        print("The unit system is the following:")
        for u in this.U:
            if u!="none":
                print("\t" + u + " is expressed in " + this.U[u])
    def make_umbrella(this, d, k):
        this.umbrella["umbrella"][0] = True
        this.umbrella["pair_distance"][0] = d
        this.umbrella["stiffness"][0] = k
    # Main function, used only from the progpoly class "pp"
    def _ProgPoly__launch(this, pp):
        ps = this.default_parameters
        this.__speak("Changing the parameter in the nvt.mdp file")
        pp.change("template/nvt.mdp", pp.folder+"nvt.mdp", ["continuation", "duration_step", "number_step", "frequency_output_x", "frequency_output_v", "frequency_output_f", "frequency_output_e", "frequency_compressed_output", "temperature", "generate_velocities", "SEED", "field", "field_z"], ["no", ps["duration_step"][0], ps["number_step"][0], ps["frequency_output_x"][0], ps["frequency_output_v"][0], ps["frequency_output_f"][0], ps["frequency_output_e"][0], ps["frequency_compressed_output"][0], ps["temperature"][0], "yes", ps["SEED"][0], ps["field"][0], "0"])
        if this.umbrella["umbrella"][0]:
            with open(pp.folder+"nvt.mdp", "a") as f:
                f.write("\npull = yes\npull_ncoords = 1\npull_ngroups = 2\npull_group1_name = IP\npull_group2_name = IM\npull_coord1_type = umbrella\npull_coord1_geometry = distance\npull_coord1_dim = Y Y N\npull_coord1_groups = 1 2\npull_coord1_k = "+str(this.umbrella["stiffness"][0])+"\npull_coord1_init = "+str(this.umbrella["pair_distance"][0])+"\npull_coord1_rate = 0\n")
        pp.call("gmx grompp -f "+pp.folder+"nvt.mdp -c "+pp.folder+"em.gro -r "+pp.folder+"em.gro -p "+pp.folder+"topol.top -n "+pp.folder+"index.ndx -o "+pp.folder+"nvt.tpr")
        if ps["number_step"][0]>0:
            if not pp.call("gmx mdrun -v -deffnm "+pp.folder+"nvt --nt " + str(pp.ncore) + " -pin on -pinoffset " + str(pp.offset) + " -pinstride 2"):
                return False
        else:
            subprocess.run("cp "+pp.folder+"em.gro "+pp.folder+"nvt.gro", shell=True)
        return True

