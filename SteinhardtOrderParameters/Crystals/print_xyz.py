import bcc

Nx = 6
Ny = 6
Nz = 6
xA = 1.0
diameter = 1.0

def printxyz(filename):
    file = open(filename, 'w')
    print (N, file = file)
    print ("Box: ", lx,ly,lz, file = file )
    for ip in range(N):
        if str( type[ip] ) == "1":
            print (type[ip], pos[ip,0], pos[ip,1], pos[ip,2], file = file)
    for ip in range(N):
        if str( type[ip] ) != "1":
            print (type[ip], pos[ip,0], pos[ip,1], pos[ip,2], file = file)
    file.close()


## print bcc disordered:

(lx,ly,lz,pos,type)=bcc.disordered(Nx, Ny, Nz, xA, diameter)
N= len(type)
filename = './bcc_disordered.xyz'
printxyz(filename)


