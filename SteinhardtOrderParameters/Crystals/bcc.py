####################################################################
######  Module to construct different types of bcc crystals
######  Hasan Zerze, April, 11th, 2017
####################################################################

## Inputs are:
        ## Nx, Ny, Nz: number of unit cells in x, y and z directions, respectively
        ## xA: Fraction of A-type particles
        ## diameter: contacting particle-particle distance
## Outputs are:
        ## Lx, Ly, Lz: box size
        ## pos[N,3]: Particle positions
        ## type[N]:  Particle types, 1 or 2

def     disordered(Nx, Ny, Nz, xA, diameter):

    import math
    import random
    import numpy as np

    lattice = (2/math.sqrt(3))*diameter

    NN=Nx*Ny*2*Nz  # particle number is 2 in unit cell
    NA = int( xA * NN )
    NB = NN - NA

    Lx = Nx * lattice
    Ly = Ny * lattice
    Lz = Nz * lattice

    pos = np.empty([NN,3], dtype=float)
    type = np.empty([NN], dtype=int)

    a1_x=0.0; a1_y=0.0; a1_z=0.0   # He
    a2_x=lattice/2; a2_y=lattice/2; a2_z=lattice/2  # H

    iatom = 0
    iA = 0
    iB = 0
    for i in range (1, Nx+1):
        for j in range (1,Ny+1):
            for k in range (1,Nz+1):
                pos[iatom,0] = a1_x + lattice*(i-1)
                pos[iatom,1] = a1_y + lattice*(j-1)
                pos[iatom,2] = a1_z + lattice*(k-1)
                xrandom = random.random()
                if (float(NA - iA)/float(NN - iA -iB)) > xrandom:
                    type[iatom] = 1
                    iA = iA + 1
                else:
                    type[iatom] = 2
                    iB = iB +1
                iatom = iatom + 1

    for i in range (1, Nx+1):
        for j in range (1,Ny+1):
            for k in range (1,Nz+1):
                pos[iatom,0] = a2_x + lattice*(i-1)
                pos[iatom,1] = a2_y + lattice*(j-1)
                pos[iatom,2] = a2_z + lattice*(k-1)
                xrandom = random.random()
                if (float(NA - iA)/float(NN - iA -iB)) > xrandom:
                    type[iatom] = 1
                    iA = iA + 1
                else:
                    type[iatom] = 2
                    iB = iB +1

                iatom = iatom + 1

    return Lx, Ly, Lz, pos, type

def     cscl(Nx, Ny, Nz, diameter):

    import math
    import random
    import numpy as np

    lattice = (2/math.sqrt(3))*diameter
    NN=Nx*Ny*2*Nz  # particle number is 2 in unit cell

    Lx = Nx * lattice
    Ly = Ny * lattice
    Lz = Nz * lattice

    pos = np.empty([NN,3], dtype=float)
    type = np.empty([NN], dtype=int)

    a1_x=0.0; a1_y=0.0; a1_z=0.0   # He
    a2_x=lattice/2; a2_y=lattice/2; a2_z=lattice/2  # H

    iatom = 0
    for i in range (1, Nx+1):
        for j in range (1,Ny+1):
            for k in range (1,Nz+1):
                pos[iatom,0] = a1_x + lattice*(i-1)
                pos[iatom,1] = a1_y + lattice*(j-1)
                pos[iatom,2] = a1_z + lattice*(k-1)
                type[iatom] = 1
                iatom = iatom + 1

    for i in range (1, Nx+1):
        for j in range (1,Ny+1):
            for k in range (1,Nz+1):
                pos[iatom,0] = a2_x + lattice*(i-1)
                pos[iatom,1] = a2_y + lattice*(j-1)
                pos[iatom,2] = a2_z + lattice*(k-1)
                type[iatom] = 2
                iatom = iatom + 1

    return Lx, Ly, Lz, pos, type

