import numpy as np
import matplotlib.pyplot as plt

# Create a hcp lattice and return it as 3 files of size Nw, Np, Nm randomly
def gen_hcp(d, xi, yi, zi, xf, yf, zf, Nw, Np, Nm, SEED):
    ypar = 1
    zpar = 1

    tab = []
    def try_write(x, y, z, tab):
        if x>=xi and y>=yi and z>=zi and x<xf and y<yf and z<zf:
            tab.append([x, y, z])
    dx = d
    dy = 3**0.5*d/2
    dz = (2/3)**0.5*d
    Nx = int((xf-xi)/dx)
    Ny = int((yf-yi)/dy)
    Nz = int((zf-zi)/dz)
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                if k%2 == 0:
                    if j%2 == 0:
                        try_write(xi+i*dx, yi+j*dy, zi+k*dz, tab)
                    else:
                        try_write(xi+i*dx+dx/2, yi+j*dy, zi+k*dz, tab)
                else:
                    if j%2 == 0:
                        try_write(xi+i*dx+dx/2, yi+j*dy+dy/2, zi+k*dz, tab)
                    else:
                        try_write(xi+i*dx, yi+j*dy+dy/2, zi+k*dz, tab)
    tab = np.array(tab)
    print(tab[:,2].min())
    np.random.seed(SEED)
    np.random.shuffle(tab)
    tab_w = tab[0:Nw]
    tab_p = tab[Nw:Nw+Np]
    tab_m = tab[Nw+Np:Nw+Np+Nm]
    return tab_w, tab_p, tab_m

# Create a hcp lattice and return it as 3 files of size Nw, Np, Nm randomly. Center the lattice on the z axis
def gen_hcp_cz(d, xi, yi, zi, xf, yf, zf, Nw, Np, Nm, SEED):
    ypar = 1
    zpar = 1

    tab = []
    def try_write(x, y, z, tab):
        if x>=xi and y>=yi and z>=zi and x<xf and y<yf and z<zf:
            tab.append([x, y, z])
    dx = d
    dy = 3**0.5*d/2
    dz = (2/3)**0.5*d
    Nx = int((xf-xi)/dx)
    Ny = int((yf-yi)/dy)
    Nz = int((zf-zi)/dz)
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                if k%2 == 0:
                    if j%2 == 0:
                        try_write(xi+i*dx, yi+j*dy, zi+k*dz, tab)
                    else:
                        try_write(xi+i*dx+dx/2, yi+j*dy, zi+k*dz, tab)
                else:
                    if j%2 == 0:
                        try_write(xi+i*dx+dx/2, yi+j*dy+dy/2, zi+k*dz, tab)
                    else:
                        try_write(xi+i*dx, yi+j*dy+dy/2, zi+k*dz, tab)
    tab = np.array(tab)
    tab[:,2] += (zf-zi)/2-Nz*dz/2
    np.random.seed(SEED)
    np.random.shuffle(tab)
    tab_w = tab[0:Nw]
    tab_p = tab[Nw:Nw+Np]
    tab_m = tab[Nw+Np:Nw+Np+Nm]
    return tab_w, tab_p, tab_m

def gen_hcp_cz_excludecircle(d, xi, yi, zi, xf, yf, zf, Nw, Np, Nm, x0, y0, R, SEED):
    tab = []
    def try_write(x, y, z, tab):
        if x>=xi and y>=yi and z>=zi and x<xf and y<yf and z<zf and (x-x0)**2+(y-y0)**2>R**2:
            tab.append([x, y, z])
    dx = d
    dy = 3**0.5*d/2
    dz = (2/3)**0.5*d
    Nx = int((xf-xi)/dx)
    Ny = int((yf-yi)/dy)
    Nz = int((zf-zi)/dz)
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                if k%2 == 0:
                    if j%2 == 0:
                        try_write(xi+i*dx, yi+j*dy, zi+k*dz, tab)
                    else:
                        try_write(xi+i*dx+dx/2, yi+j*dy, zi+k*dz, tab)
                else:
                    if j%2 == 0:
                        try_write(xi+i*dx+dx/2, yi+j*dy+dy/2, zi+k*dz, tab)
                    else:
                        try_write(xi+i*dx, yi+j*dy+dy/2, zi+k*dz, tab)
    tab = np.array(tab)
    tab[:,2] += (zf-zi)/2-Nz*dz/2
    np.random.seed(SEED)
    np.random.shuffle(tab)
    tab_w = tab[0:Nw]
    tab_p = tab[Nw:Nw+Np]
    tab_m = tab[Nw+Np:Nw+Np+Nm]
    return tab_w, tab_p, tab_m