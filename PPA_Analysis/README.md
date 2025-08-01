
################################
## Compilation:
################################

Make sure you have the XDRFILE library installed (used to read `.xtc` files).
Then, to compile the code, enter the following command:

g++ -O3 main.cpp config.cpp xtc_reader.cpp ppa_solver.cpp -o run_ppa -I./libxdrfile/include -L./libxdrfile/build/lib -lxdrfile

