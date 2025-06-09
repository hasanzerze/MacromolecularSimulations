# Bead-Spring Polymer Chain Generator for LAMMPS

## Overview

This Python script generates LAMMPS-compatible data files for systems composed of bead-spring polymer chains. The chains are constructed as random walks in 3D, and the data file includes atoms, bonds, and molecule information formatted according to LAMMPS standards.

**Author:** Hasan Zerze  
**Last Updated:** May 16th, 2025

---

## Features

- Generates linear bead-spring chains with a fixed bond length.
- Automatically places each chain randomly in a cubic simulation box.
- Writes a `.data` file for LAMMPS containing:
  - Atoms
  - Bonds
  - Masses
  - Box size

---

## Prerequisites

- Python 3.x
- `numpy`

To install required dependencies:

```bash
pip install numpy

## Input: `system_setup.dat`
This file contains the system setup parameters in plain text format. Each line should define one variable, for example:

Nchain 100
Natom_per_chain 10
bondlength 1.0
Lbox 50.0

Where:

Nchain: Number of polymer chains

Natom_per_chain: Number of atoms per chain

bondlength: Fixed bond length between adjacent beads

Lbox: Size of the cubic simulation box (x=Lbox, y=Lbox, z=Lbox)

## Output: `file.data`
The LAMMPS-compatible data file containing:

Number of atoms and bonds
Atom types, bond types
Box dimensions
Masses
Atom coordinates
Bond connectivity
Each atom is tagged with a molecule ID and atom type (default: 1). Bonds are added between consecutive atoms in each chain.

## Usage
Run the Python script in the same directory as your system_setup.dat file:

python main_setup.py

The output will be saved as file.data in the current directory.

## Notes
- Chains are generated as random walks by random sampling directions on a unit sphere and ensuring they stay within the bond length constraint.
- No angle or dihedral interactions are included in this script, making it suitable for coarse-grained, freely-jointed chain simulations.

## License
This project is released under the MIT License.

## Citation
If you use this code for your work or publication, please cite as:

@misc{zerze_beadspring_2025,
  author       = {Hasan Zerze},
  title        = {Bead-Spring Polymer Chain Generator for LAMMPS},
  year         = {2025},
  howpublished = {\url{https://github.com/hasanzerze/MacromolecularSimulations}},
  note         = {Python script for generating linear bead-spring chains}
}


