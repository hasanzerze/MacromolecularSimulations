#pragma once
#include <vector>

struct Position {
    int x, y;
    Position(int x_=0, int y_=0) : x(x_), y(y_) {}
    bool operator==(const Position& other) const {
        return x == other.x && y == other.y;
    }
};

class Lattice {
public:
    Lattice(int size);
    bool isOccupied(int x, int y) const;
    void occupy(int x, int y);
    void clear();

private:
    int size;
    std::vector<std::vector<bool>> grid;
};

