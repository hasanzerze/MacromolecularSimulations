#include "lattice.h"

Lattice::Lattice(int size_) : size(size_), grid(size_, std::vector<bool>(size_, false)) {}

bool Lattice::isOccupied(int x, int y) const {
    return grid[x][y];
}

void Lattice::occupy(int x, int y) {
    grid[x][y] = true;
}

void Lattice::clear() {
    for (auto& row : grid) std::fill(row.begin(), row.end(), false);
}
