import numpy as np

Npar = 500
kB = 1.0
Temp = 0.85

data = np.loadtxt("thermo_output.txt", dtype=float, comments='#', usecols=(1))

Ndata = np.shape(data)[0]
Nstart = int(Ndata / 2)

print (Ndata, Nstart)
energy = data[Nstart:]
Emean = np.mean(energy)
E2mean = np.mean(energy * energy)

Cv = 3.0 * kB / 2.0 + ( E2mean - Emean * Emean ) / Npar / kB / Temp / Temp

print (Ndata, Emean, Cv)
#print (data)


