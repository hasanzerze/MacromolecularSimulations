#!/bin/bash

tail --lines=+203 lammps.out | head --lines=1003 > thermo_nvt.dat

tail --lines=+1242 lammps.out | head --lines=1002 > thermo_npt.dat

