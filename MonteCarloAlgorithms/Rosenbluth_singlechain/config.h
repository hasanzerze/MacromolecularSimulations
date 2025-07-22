#pragma once
#include <string>

namespace config {
    extern int LATTICE_SIZE;
    extern int CHAIN_LENGTH;
    extern double EPSILON;
    extern int NUM_TRIALS;
    extern int SEED;

    void loadConfigFromFile(const std::string& filename);
}

