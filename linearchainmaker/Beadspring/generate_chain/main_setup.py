##  Hasan Zerze, Last updated: May 16th, 2025
##  The code to generate LAMMPS data file for
##      bead-spring chains

import numpy as np
import math
import random
from numpy import random

def make_chain(bondlength, Natom_per_chain):
    pos = np.zeros((Natom_per_chain, 3), dtype=float)
    for i in range(1,Natom_per_chain):
        direction = [1., 1., 1.]
        while ( np.linalg.norm( direction ) > 1.0 ):
            direction = 2. * random.random(size=3) - 1.
        
        segment = direction / np.linalg.norm( direction )
        segment *= bondlength

        pos[i, :] = pos[i-1, :] + segment            

    type_seq = Natom_per_chain * [1]

    return pos, type_seq

def translate_molecule(pos, shift):
    Natom_per_mol = np.shape(pos)[0]
    dim = np.shape(pos)[1]
    translated_pos = np.zeros( (Natom_per_mol, dim), \
                    dtype=float)

    for i in range ( Natom_per_mol ):
        translated_pos[i,:] = pos[i,:] + shift[:]

    return translated_pos

setup=open("system_setup.dat","r")

for line in setup.readlines():
        data = line.split()
        vars()[data[0]] = float(data[1])

### Estimate the total number of atoms:
Nchain = int(Nchain)
Natom_per_chain = int(Natom_per_chain)
Natoms = Nchain * Natom_per_chain

### Estimate all the positions:
positions = np.zeros((Natoms,3), dtype=float)
moltag = np.zeros((Natoms), dtype=int)
atomtype = np.zeros((Natoms), dtype=int)
iatom = 0
imol = 1

Nbonds_per_chain = Natom_per_chain - 1
Nbonds = Nchain * Nbonds_per_chain

if Nchain != 0:
    for i in range(Nchain):
        (Chain_ref_pos, type_seq_Chain) \
            = make_chain(bondlength, Natom_per_chain)

        (istart, istop) = (iatom, iatom + Natom_per_chain)
        randpos = random.random(size=3) * Lbox
        positions[istart:istop,:] = translate_molecule( \
                        Chain_ref_pos, randpos)

        atomtype[istart:istop] = type_seq_Chain
        moltag[istart:istop] = Natom_per_chain * [imol]
        iatom += Natom_per_chain
        imol += 1

## Producing the data file:
Natomtypes = 1
Nbondtypes = 1

mass = [1.0]

datafile = open("file.data", "w")

print ("LAMMPS data file for linear chains", file=datafile)
print ("", file=datafile )
print (Natoms, "  atoms", file=datafile )
print (Nbonds, "  bonds", "#", file=datafile )
print (0, "  angles", file=datafile )
print (0, "  dihedrals", file=datafile )
print (0, "  impropers", file=datafile )
print ("", file=datafile )
print (Natomtypes, "  atom types", file=datafile )
print (Nbondtypes, "  bond types", file=datafile )
print (0, "  angle types", file=datafile )
print (0, "  dihedral types", file=datafile )
print (0, "  improper types", file=datafile )
print ("", file=datafile )
print (0, Lbox, "  xlo xhi", file=datafile )
print (0, Lbox, "  ylo yhi", file=datafile )
print (0, Lbox, "  zlo zhi", file=datafile )
print ("", file=datafile )
print ("Masses", file=datafile )
print ("", file=datafile )
for i in range (Natomtypes):
    print (i+1, mass[i], file=datafile )

print ("", file=datafile )
print ("Atoms", file=datafile )
print ("", file=datafile )
for i in range(Natoms):
    print (i+1, moltag[i], atomtype[i], \
            positions[i,0], positions[i,1], \
            positions[i,2], file=datafile )

print ("", file=datafile )
print ("Bonds", file=datafile )
print ("", file=datafile )
ib = 1

if Nchain != 0:
    for i in range( Nchain ):
        # calculate the index of the starting atom
        k = i * Natom_per_chain

        for j in range( Nbonds_per_chain ):
            k += 1
            print ( ib, 1, k, k+1, file=datafile )  
            ib += 1


