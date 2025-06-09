# Monte Carlo Simulation of Lennard-Jones Particles

This project implements a canonical ensemble (NVT) **Monte Carlo simulation** of particles interacting via the **Lennard-Jones (LJ) potential** in 3D space using **Metropolis sampling** and **cell list optimization** for neighbor searching.

## Features

- Simulation of monoatomic particles in a cubic periodic box
- Uses **Lennard-Jones interactions** with cutoff and long-range corrections
- **Periodic boundary conditions** and **minimum image convention**
- **Metropolis algorithm** for sampling particle configurations
- Computes:
  - Total energy
  - Pressure (via virial expression and long-range correction)
  - Acceptance ratio
- Implements **linked cell list** to improve performance of neighbor searches

---

## Getting Started

### Prerequisites

- C++ compiler supporting C++11 or higher
- Standard C++ libraries (no external dependencies)

### Building the Code

Compile using any C++ compiler. For example:

```bash
g++ -std=c++11 -O2 -o lj_mc_sim main.cpp

```
## Running the Simulation:

```bash
./lj_mc_sim
```

## Parameters:

```
Parameter	Description	Value
NUM_PARTICLES	Number of particles	500
BOX_LENGTH	Length of cubic simulation box	8.22
RCUT	Cutoff radius for LJ interaction	3.0
STEP_SIZE	Maximum displacement in a trial move	0.1
temperature	Simulation temperature (in LJ units)	0.85
LJ_EPSILON	LJ energy scale	1.0
LJ_SIGMA	LJ length scale	1.0
```
## Code Structure:

main() initializes the simulation and runs it.
Simulation class contains all logic for particle motion, energy calculation, and pressure estimation.
Particle stores the position and current cell index of each particle.
Vector3d is a helper structure for 3D vector operations.

## Notes:

The initial particle configuration avoids overlaps using a rejection method.
Long-range energy and pressure corrections are applied using analytical integrals.
The simulation assumes reduced LJ units: ε = σ = kB = 1.

## License:
This code is provided for educational and research purposes. Feel free to modify and extend it.

## Author
Hasan Zerze

