# External library
import numpy as np
import subprocess

# Internal library
import initialisation
import thermalisation
import simulation
import progpoly

def pin(i):
    Nslots=4
    nt=4
    ntotalcores = Nslots*nt
    return np.r_[np.arange(0,ntotalcores,2*nt),
        np.arange(0,ntotalcores,2*nt)+1][:Nslots].astype(int)[i]
'''
# Code for transport simulation
if __name__ == "__main__":
    ## Initialization parameters
    SEED = 1
    for E in 10.**(np.array([-1])):
        h=0.7
        a = 0.142
        dx = 3**0.5*a
        dy = 3*a
        Lx = 20
        Ly = 20
        Lx_corr = int(Lx/dx)*dx
        Ly_corr = int(Ly/dy)*dy
        dens = 11.8 #11.2 # Default 11.2 # 2D:21.22
        Nw = int(dens*Lx_corr*Ly_corr)
        print(Nw)
        init = initialisation.Initialisation()
        init.fix(charge_p=0, charge_m=0, length_x=Lx, length_y=Ly, length_z=20, number_w=Nw, number_p=0, number_m=0, SEED=SEED, height=h)
        ## Thermalization parameters
        therm = thermalisation.Thermalisation()
        therm.fix(number_step = 1000000, SEED=SEED, field=E)
        ## Simulation parameters
        run = simulation.Simulation()
        run.fix(frequency_output_x=10000, frequency_output_v=0, frequency_output_f=0, frequency_compressed_output=0, frequency_output_e=0, field=E, number_step = 50000000, SEED=SEED)
        
        ## Simulation
        pp = progpoly.ProgPoly("/usr/local/gromacs/bin/gmx", init, therm, run, 16, 0, True, True, True)
        pp.name(SEED=SEED, h=h, dens=dens, E=E, time="50ns")
        pp.launch()
    #pp.clean()
    #pp.compute_cluster()
    #pp.call("gmx trjconv -f "+pp.folder+"run.trr -s "+pp.folder+"run.tpr -n "+pp.folder+"index.ndx -o simu.pdb -skip 100")
    ##pp.call("/usr/local/gromacs/bin/gmx trjconv -f run.trr -s run.tpr -n index.ndx -o simu.pdb -skip 100")
    #pp.clean()
'''
# Code for transport simulation
if __name__ == "__main__":
    ## Initialization parameters
    SEED = 1
    h=0.7
    a = 0.142
    dx = 3**0.5*a
    dy = 3*a
    Lx = 20
    Ly = 20
    Lx_corr = int(Lx/dx)*dx
    Ly_corr = int(Ly/dy)*dy
    dens = 11.8 #11.2 # Default 11.2 # 2D:21.22
    Nw = int(dens*Lx_corr*Ly_corr)
    print(Nw)
    init = initialisation.Initialisation()
    init.fix(charge_p=0, charge_m=0, length_x=Lx, length_y=Ly, length_z=20, number_w=Nw, number_p=0, number_m=0, SEED=SEED, height=h)
    ## Thermalization parameters
    therm = thermalisation.Thermalisation()
    therm.fix(number_step = 1, SEED=SEED, field=0)
    ## Simulation parameters
    run = simulation.Simulation()
    run.fix(frequency_output_x=10000, frequency_output_v=0, frequency_output_f=0, frequency_compressed_output=0, frequency_output_e=0, field=0, number_step = 51000000, SEED=SEED)
    
    ## Simulation
    pp = progpoly.ProgPoly("/usr/local/gromacs/bin/gmx", init, therm, run, 8, 0, True, True, True)
    pp.name(SEED=SEED, h=h, dens=dens, time="50ns_nothermal")
    pp.launch()
    #pp.clean()
    #pp.compute_cluster()
    #pp.call("gmx trjconv -f "+pp.folder+"run.trr -s "+pp.folder+"run.tpr -n "+pp.folder+"index.ndx -o simu.pdb -skip 100")
    ##pp.call("/usr/local/gromacs/bin/gmx trjconv -f run.trr -s run.tpr -n index.ndx -o simu.pdb -skip 100")
    #pp.clean()
