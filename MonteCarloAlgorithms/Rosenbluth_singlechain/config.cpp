#include "config.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>

namespace config {
    int LATTICE_SIZE = 0;
    int CHAIN_LENGTH = 0;
    double EPSILON = 0.0;
    int NUM_TRIALS = 0;
    int SEED = 0;

    void loadConfigFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Error: Could not open config file: " << filename << std::endl;
            std::exit(1);
        }

        std::string key;
        double value;
        std::map<std::string, double> constants;

        while (file >> key >> value) {
            constants[key] = value;
        }

        if (constants.count("LATTICE_SIZE"))  LATTICE_SIZE = static_cast<int>(constants["LATTICE_SIZE"]);
        if (constants.count("CHAIN_LENGTH"))  CHAIN_LENGTH = static_cast<int>(constants["CHAIN_LENGTH"]);
        if (constants.count("EPSILON"))       EPSILON = constants["EPSILON"];
	if (constants.count("NUM_TRIALS"))       NUM_TRIALS = constants["NUM_TRIALS"];
	if (constants.count("SEED"))       SEED = constants["SEED"];
    }
}

