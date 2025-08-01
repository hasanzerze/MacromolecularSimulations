#include "config.h"
#include <fstream>
#include <sstream>
#include <iostream>

namespace config {
    std::string XTC_FILE;
    int NUM_ATOMS = 0;
    int FRAME_STRIDE = 1;

    int NUM_CHAINS = 0;
    int MONOMERS_PER_CHAIN = 0;

    std::vector<Bond> BONDS;
    std::vector<std::vector<int>> CHAINS;

    LJParams LJ_PARAMS;
    PPASolverParams PPA_PARAMS;

    void loadConfigFromFile(const std::string& filename) {
        std::ifstream in(filename);
        if (!in) {
            std::cerr << "Error opening config file: " << filename << std::endl;
            exit(1);
        }

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            std::string key;
            iss >> key;

            if (key == "XTC_FILE") {
                iss >> XTC_FILE;
            } else if (key == "NUM_ATOMS") {
                iss >> NUM_ATOMS;
            } else if (key == "FRAME_STRIDE") {
                iss >> FRAME_STRIDE;
            } else if (key == "LJ_EPSILON") {
                iss >> LJ_PARAMS.epsilon;
            } else if (key == "LJ_SIGMA") {
                iss >> LJ_PARAMS.sigma;
            } else if (key == "PPA_MAX_ITER") {
                iss >> PPA_PARAMS.max_iterations;
            } else if (key == "PPA_BOND_TOL") {
                iss >> PPA_PARAMS.bond_length_tolerance;
            } else if (key == "PPA_TIMESTEP") {
                iss >> PPA_PARAMS.timestep;
	    } else if (key == "PPA_BOND_CONST") {
                iss >> PPA_PARAMS.ppa_bond_const;
            } else if (key == "NUM_CHAINS") {
	        iss >> NUM_CHAINS;
//	    } else if (key == "DEBUG_OUTPUT:") {
//	        int val;
//	        iss >> val;
//	        PPA_PARAMS.debug_output = (val != 0);
	    } else if (key == "MONOMERS_PER_CHAIN") {
	        iss >> MONOMERS_PER_CHAIN;
	    } else {
                std::cerr << "Unknown config key: " << key << std::endl;
            }
        }

        if (NUM_ATOMS != NUM_CHAINS * MONOMERS_PER_CHAIN) {
            std::cerr << "Error: NUM_ATOMS must equal NUM_CHAINS × MONOMERS_PER_CHAIN.\n";
            std::exit(1);
        }

	config::CHAINS.clear();
	for (int c = 0; c < config::NUM_CHAINS; ++c) {
	    std::vector<int> chain;
	    for (int m = 0; m < config::MONOMERS_PER_CHAIN; ++m) {
	        chain.push_back(c * config::MONOMERS_PER_CHAIN + m);
	    }
	    config::CHAINS.push_back(chain);
	}

        // Now generate bonds automatically
        generateBonds();
    }
    void generateBonds() {
        BONDS.clear();  // start fresh

        for (int c = 0; c < NUM_CHAINS; ++c) {
            int offset = c * MONOMERS_PER_CHAIN;
            for (int i = 0; i < MONOMERS_PER_CHAIN - 1; ++i) {
                BONDS.push_back({offset + i, offset + i + 1});
            }
        }
	std::cout << "[INFO] Bonds generated: " << BONDS.size() << std::endl;
    }
}

