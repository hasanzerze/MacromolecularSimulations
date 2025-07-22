#include <random>
#include <iostream>
#include "config.h"
#include "lattice.h"
#include "chain.h"

int main() {
    config::loadConfigFromFile("input.dat");

    std::cout << "Simulation Parameters:\n";
    std::cout << "Lattice Size: " << config::LATTICE_SIZE << "\n";
    std::cout << "Chain Length: " << config::CHAIN_LENGTH << "\n";
    std::cout << "Epsilon: " << config::EPSILON << "\n";
    std::cout << "Number of Trials: " << config::NUM_TRIALS << std::endl;
    std::cout << "Random Seed: " << config::SEED << std::endl;

    std::mt19937 rng(config::SEED);

    Lattice lattice(config::LATTICE_SIZE);

    double totalWeight = 0.0;
    double totalWeightSquared = 0.0;
    double totalContacts = 0.0;
    int success = 0;
    double totalRee2 = 0.0;
    double totalRg2 = 0.0;

    for (int trial = 0; trial < config::NUM_TRIALS; ++trial) {
        Chain chain(config::CHAIN_LENGTH);
        double weight = chain.grow(lattice, config::EPSILON, rng);
        if (weight > 0.0) {
            totalWeight += weight;
            success++;
	    totalWeightSquared += weight * weight;

	    // Count total contacts in final chain
	    double contacts = 0;
	    const auto& positions = chain.getPositions();

	    for (size_t i = 0; i < positions.size(); ++i) {
	        contacts += chain.countContacts(i, lattice);
	    }
            // Each contact counted twice
            totalContacts += contacts / 2.0;

	    // End-to-end distance squared
            const Position& p0 = positions.front();
            const Position& pN = positions.back();
            double dx = pN.x - p0.x;
            double dy = pN.y - p0.y;
            double ree2 = dx * dx + dy * dy;
            totalRee2 += ree2;

            // Radius of gyration squared
            double xcm = 0.0, ycm = 0.0;
            for (const auto& pos : positions) {
                xcm += pos.x;
                ycm += pos.y;
            }
            xcm /= positions.size();
            ycm /= positions.size();

            double rg2 = 0.0;
            for (const auto& pos : positions) {
                double dx = pos.x - xcm;
                double dy = pos.y - ycm;
                rg2 += dx * dx + dy * dy;
            }
            rg2 /= positions.size();
            totalRg2 += rg2;

	} //else {
        // Print configuration of failed chain
	//    std::cout << "Trial " << trial << " failed to grow full chain.\n";
        //    const auto& positions = chain.getPositions();
        //    std::cout << "Chain length reached: " << positions.size() << "\n";
        //    for (size_t i = 0; i < positions.size(); ++i) {
        //        std::cout << "  Bead " << i << ": (" << positions[i].x << ", " << positions[i].y << ")\n";
        //    }
        //    std::cout << "-----\n";
        //}
    }

    if (success > 0) {
        double avgWeight = totalWeight / success;
        double varWeight = (totalWeightSquared / success) - (avgWeight * avgWeight);
        double stdWeight = std::sqrt(varWeight);
        double avgContacts = totalContacts / success;

        std::cout << "Successful chains: " << success << " / " << config::NUM_TRIALS << "\n";
        std::cout << "Average Rosenbluth weight: " << avgWeight << " ± " << stdWeight << "\n";
        std::cout << "Average number of contacts: " << avgContacts << "\n";
	std::cout << "Average end-to-end distance squared: " << totalRee2 / success << "\n";
	std::cout << "Average radius of gyration squared: " << totalRg2 / success << "\n";

    } else {
        std::cout << "All trials failed. Try increasing lattice size or decreasing chain length.\n";
    }

    return 0;
}


