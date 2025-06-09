# Steinhardt Order Parameter Calculator

This C++ program computes the **Steinhardt bond orientational order parameters** (specifically \( q_l \)) from an XTC trajectory file. These parameters are commonly used to quantify local structural order in molecular or particulate systems.

---

## Project Structure

- `steinhardt.cpp`: Main program to compute \( q_l \) from an XTC trajectory.
- `libxdrfile/`: Directory containing the locally built `xdrfile` library.
- `README.md`: Documentation file.

---

## Requirements

- A C++ compiler supporting C++17 (e.g., `g++`)
- The `xdrfile` library built and located in:
  - Headers: `libxdrfile/include`
  - Compiled library: `libxdrfile/build/lib/libxdrfile.a` or `.so`

---

## Compilation

Use the following command to compile:

g++ -o steinhardt steinhardt.cpp -I./libxdrfile/include -L./libxdrfile/build/lib -lxdrfile

## Notes:
Make sure you build the xdrfile library first. If using CMake:

cd libxdrfile
mkdir -p build && cd build
cmake ..
make

## Usage:

./steinhardt <xtc file>

Example xtc files for some crystal structures are available under the directory named `Crystals`

## Output
Prints the scalar order parameter ql for each atom in each frame.

Example output:

Time: 0 ps Box size = 20, 20, 20
ql = 0.4893
ql = 0.5121
...

## Parameters
These can be modified directly in the source code:
l = 4: Order of the spherical harmonics.
r_neigh = 1.4: Cutoff radius for neighbor search (in reduced units).

## References
P. J. Steinhardt, D. R. Nelson, M. Ronchetti, Phys. Rev. B, 28, 784 (1983).

## Author
Hasan Zerze


