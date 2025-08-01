#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
#include "vector3d.h"

namespace config {
    struct Bond { int i, j; };

    struct LJParams {
        double epsilon;
        double sigma;
    };

    struct PPASolverParams {
        double timestep;
        double bond_length_tolerance;
	double ppa_bond_const;
        int max_iterations;
        bool debug_output = false;
    };

    extern std::string XTC_FILE;
    extern int NUM_ATOMS;
    extern int NUM_CHAINS;
    extern int MONOMERS_PER_CHAIN;
    extern int FRAME_STRIDE;
    extern std::vector<Bond> BONDS;
    extern std::vector<std::vector<int>> CHAINS;
    extern LJParams LJ_PARAMS;
    extern PPASolverParams PPA_PARAMS;

    void loadConfigFromFile(const std::string& filename);
    void generateBonds();
}

#endif

