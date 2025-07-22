#pragma once
#include "lattice.h"
#include <vector>
#include <random>

class Chain {
public:
    Chain(int length);
    double grow(Lattice& lattice, double epsilon, std::mt19937& rng);

    const std::vector<Position>& getPositions() const;
    int countContacts(int idx, const Lattice& lattice) const;
    int countContacts(const Position& pos, const Lattice& lattice) const;

private:
    int length;
    std::vector<Position> positions;
    double weight;

    std::vector<Position> getFreeNeighbors(const Position& pos, const Lattice& lattice) const;
};
