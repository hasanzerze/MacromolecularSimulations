################################
# Hasan Zerze, 08/01/2025
# Primitive Path Analysis (PPA) for Polymer Entanglement
################################

This repository contains a C++ implementation of the Primitive Path Analysis (PPA)
algorithm used to study topological entanglements in polymer systems.
It reads atomic trajectories from GROMACS xtc files and performs energy minimization
to compute quantities such as entanglement length, primitive path length, and radius of gyration.

################################
## Purpose
################################

The goal of this project is to:
- Analyze entanglement properties of polymer chains.
- Compute statistical observables like primitive path length (Lpp), 
  mean squared end-to-end distance (Ree2), and entanglement length (Ne).
- Provide a reference implementation for energy minimization using 
  the FIRE algorithm.

################################
## Features
################################

- Reads xtc files using the xdrfile library.
- Supports the use of custom simulation parameters via config files.
- Performs iterative minimization per frame with desired frame read frequency.
- Allows multiple rounds of minimization to monitor convergence of entanglement statistics.
- Outputs minimized coordinates back into an xtc formatted file.

################################
## Dependencies
################################

- C++17 or newer
- xdrfile library (linked during compilation)

################################
## Compilation:
################################

Make sure you have the XDRFILE library installed (used to read xtc files).
Then, to compile the code, enter the following command:

g++ -O3 main.cpp config.cpp xtc_reader.cpp ppa_solver.cpp -o run_ppa -I./libxdrfile/include -L./libxdrfile/build/lib -lxdrfile

################################
## Input file:
################################
The input file should be named as "input.dat"
An example content for the "input.dat" is copied below:

XTC_FILE traj.xtc
FRAME_STRIDE 1000

NUM_ATOMS 20000
NUM_CHAINS 200
MONOMERS_PER_CHAIN 100

# LJ parameters
LJ_EPSILON 0.1
LJ_SIGMA 0.4

# PPA solver parameters
PPA_MAX_ITER 500
PPA_BOND_TOL 1e-4
PPA_BOND_CONST 30.0
PPA_TIMESTEP 2e-4


################################
## How to run:
################################

Once the input.dat file is prepared as above, run using the following command:

./run_ppa

################################
## Observables computed
################################

Primitive Path Length, Lpp
Mean Squared End-to-End Distance, Ree2
Kuhn length of the primitive chain, app = Ree2 / Lpp
Bond distance of the primitive chain, bpp = Lpp / natoms_per_chain
Entanglement Length, Ne = app / bpp

################################
## References
################################

Everaers, R., Sukumaran, S. K., Grest, G. S., Svaneborg, C., Sivasubramanian, A., & Kremer, K. (2004). Rheology and microscopic topology of entangled polymeric liquids. Science, 303(5659), 823–826.

Bitzek, E., Koskinen, P., Gähler, F., Moseler, M., & Gumbsch, P. (2006). Structural Relaxation Made Simple. Physical Review Letters, 97(17), 170201. (for FIRE algorithm)

################################


