import bcc

Nx = 6
Ny = 6
Nz = 6
xA = 1.0
diameter = 1.0

def printdatafile(filename):
    file = open(filename, 'w')

    print ("LAMMPS Description           (1st line of file)", file = file)
    print (" ", file = file)
    print (N, "atoms         (this must be the 3rd line, 1st 2 lines are ignored)", file = file)
    print("0 bonds                (# of bonds to be simulated)", file = file)
    print("0 angles               (include these lines even if number = 0)", file = file)
    print("0 dihedrals", file = file)
    print("0 impropers", file = file)
    print (" ", file = file)
    print("1 atom types           (# of nonbond atom types)", file = file)
    print("0 bond types          (# of bond types = sets of bond coefficients)", file = file)
    print("0 angle types", file = file)
    print("0 dihedral types      (do not include a bond,angle,dihedral,improper type)", file = file)
    print("0 improper types             line if number of bonds,angles,etc is 0)", file = file)
    print (" ", file = file)
    print(0.0, lx, "xlo xhi       (for periodic systems this is box size)", file = file)
    print(0.0, ly, "ylo yhi        for non-periodic it is min/max extent of atoms)", file = file)
    print(0.0, lz,  "zlo zhi       (do not include this line for 2-d simulations)", file = file)
    print (" ", file = file)
    print("Masses", file = file)
    print (" ", file = file)
    print (1, 1, file = file)
    print (" ", file = file)
    print ("Atoms", file = file)
    print (" ", file = file)
    iatom = 1
    imol = 1
    for ip in range(N):
        if str( type[ip] ) == "1":
            print (iatom, imol, type[ip], 0, pos[ip,0], pos[ip,1], pos[ip,2], file = file)
            iatom += 1
            imol += 1
    for ip in range(N):
        if str( type[ip] ) != "1":
            print (iatom, imol, type[ip], 0, pos[ip,0], pos[ip,1], pos[ip,2], file = file)
            iatom += 1
            imol += 1
    file.close()


## print bcc disordered:

(lx,ly,lz,pos,type)=bcc.disordered(Nx, Ny, Nz, xA, diameter)
N= len(type)
filename = './bcc_disordered.data'
printdatafile(filename)


