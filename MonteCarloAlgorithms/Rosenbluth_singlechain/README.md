# Hasan Zerze, Ph.D. July, 22nd, 2025
# Self-Avoiding Polymer Chain Simulation (2D Lattice, Rosenbluth Sampling)

This project implements a 2D lattice-based simulation of linear self-avoiding polymer chains using the **Rosenbluth sampling method**.

## Simulation Description

- Chains are grown bead-by-bead on a 2D square lattice.
- **Self-avoiding walk**: Beads do not overlap.
- **Bead–bead interactions** (contact energy `epsilon`) are included for non-bonded adjacent neighbors.
- **Rosenbluth weights** are used to handle sampling bias during chain growth.
- The simulation averages physical observables such as:
  - Rosenbluth weight
  - Number of non-bonded contacts
  - End-to-end distance
  - Radius of gyration

## Parameters (in `input.dat`)
LATTICE_SIZE 100
CHAIN_LENGTH 50
EPSILON -1.0
NUM_TRIALS 1000
SEED 42

## Files Structure
main.cpp # Main simulation loop and observable tracking
chain.h/.cpp # Chain class and growth logic
lattice.h/.cpp # 2D lattice representation
config.h/.cpp # Reads simulation parameters from input.dat
input.dat # Input file with simulation parameters
README.md # This file

## Notes

- Chains that fail to complete the desired length are discarded.
- Each contact is counted once (not double-counted).
- Observables are averaged only over successful chains.

## Compilation

g++ -std=c++17 main.cpp config.cpp lattice.cpp chain.cpp -o rosenbluth

## Author

Hasan Zerze (with substantial help from ChatGPT)

