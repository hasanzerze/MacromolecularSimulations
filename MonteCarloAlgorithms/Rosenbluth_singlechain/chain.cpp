#include "config.h"
#include "chain.h"
#include <random>
#include <cmath>
#include <algorithm>

Chain::Chain(int length_) : length(length_), weight(1.0) {}

std::vector<Position> Chain::getFreeNeighbors(const Position& pos, const Lattice& lattice) const {
    std::vector<Position> options;
    int dx[4] = {1, -1, 0, 0};
    int dy[4] = {0, 0, 1, -1};
    for (int i = 0; i < 4; ++i) {
        int nx = pos.x + dx[i];
        int ny = pos.y + dy[i];
        if (nx >= 0 && ny >= 0 && nx < config::LATTICE_SIZE && ny < config::LATTICE_SIZE &&
            !lattice.isOccupied(nx, ny)) {
            options.emplace_back(nx, ny);
        }
    }
    return options;
}

int Chain::countContacts(const Position& pos, const Lattice& lattice) const {
    int dx[4] = {1, -1, 0, 0};
    int dy[4] = {0, 0, 1, -1};
    int contacts = 0;

    for (int i = 0; i < 4; ++i) {
        int nx = pos.x + dx[i];
        int ny = pos.y + dy[i];

        if (nx >= 0 && ny >= 0 && nx < config::LATTICE_SIZE && ny < config::LATTICE_SIZE &&
            lattice.isOccupied(nx, ny)) {
            contacts++;
        }
    }

    return contacts;
}

int Chain::countContacts(int idx, const Lattice& lattice) const {
    const Position& pos = positions[idx];

    int dx[4] = {1, -1, 0, 0};
    int dy[4] = {0, 0, 1, -1};
    int contacts = 0;
    for (int i = 0; i < 4; ++i) {
        int nx = pos.x + dx[i];
        int ny = pos.y + dy[i];
        if (nx >= 0 && ny >= 0 && nx < config::LATTICE_SIZE && ny < config::LATTICE_SIZE &&
            lattice.isOccupied(nx, ny)) {

            // Check if it's a bonded neighbor (i-1 or i+1)
	    bool isBondedNeighbor = false;
            if (idx > 0 && positions[idx - 1].x == nx && positions[idx - 1].y == ny)
                isBondedNeighbor = true;
            if (idx < positions.size() - 1 &&
                positions[idx + 1].x == nx && positions[idx + 1].y == ny)
                isBondedNeighbor = true;

            if (!isBondedNeighbor)
                contacts++;
        }
    }
    return contacts;
}

double Chain::grow(Lattice& lattice, double epsilon, std::mt19937& rng) {
    lattice.clear();
    positions.clear();

    // Start from center
    Position current(config::LATTICE_SIZE / 2, config::LATTICE_SIZE / 2);
    positions.push_back(current);
    lattice.occupy(current.x, current.y);

    weight = 1.0;

    for (int step = 1; step < length; ++step) {
        auto options = getFreeNeighbors(current, lattice);
        if (options.empty()) return 0.0;  // Growth failed

        std::vector<double> probabilities;
        double totalWeight = 0.0;

        for (const auto& opt : options) {
            int contacts = countContacts(opt, lattice);
            double w = std::exp(-epsilon * contacts);
            probabilities.push_back(w);
            totalWeight += w;
        }

        std::uniform_real_distribution<> dist(0.0, totalWeight);
        double r = dist(rng);
        double cumulative = 0.0;
        size_t chosen = 0;
        for (; chosen < probabilities.size(); ++chosen) {
            cumulative += probabilities[chosen];
            if (r <= cumulative) break;
        }

        current = options[chosen];
        positions.push_back(current);
        lattice.occupy(current.x, current.y);
        weight *= totalWeight;
    }

    return weight;
}

const std::vector<Position>& Chain::getPositions() const {
    return positions;
}


