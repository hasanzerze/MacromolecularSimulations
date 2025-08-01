#ifndef CELL_LIST_H
#define CELL_LIST_H

#include "vector3d.h"
#include <vector>
#include <cmath>
#include <unordered_map>

// A 3D key for a grid cell
struct CellIndex {
    int x, y, z;
    bool operator==(const CellIndex& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

namespace std {
    template <>
    struct hash<CellIndex> {
        size_t operator()(const CellIndex& c) const {
            return (c.x * 73856093) ^ (c.y * 19349663) ^ (c.z * 83492791);
        }
    };
}

struct CellList {
    double cell_size;
    Vector3D box_size;
    std::unordered_map<CellIndex, std::vector<int>> cells;

    CellList(double cutoff, const Vector3D& box) : cell_size(cutoff), box_size(box) {}

    CellIndex positionToIndex(const Vector3D& pos, const Vector3D& box_size) const {
	Vector3D wrapped = wrap_position(pos, box_size);
        return {
            static_cast<int>(std::floor(wrapped.x / cell_size)),
            static_cast<int>(std::floor(wrapped.y / cell_size)),
            static_cast<int>(std::floor(wrapped.z / cell_size))
        };
    }

    void build(const std::vector<Vector3D>& positions, const Vector3D& box_size) {
        cells.clear();
        for (int i = 0; i < positions.size(); ++i) {
            CellIndex idx = positionToIndex(positions[i], box_size);
            cells[idx].push_back(i);
        }
    }

    std::vector<CellIndex> getNeighboringCells(const CellIndex& idx) const {
        std::vector<CellIndex> neighbors;
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    neighbors.push_back({idx.x + dx, idx.y + dy, idx.z + dz});
                }
            }
        }
        return neighbors;
    }
};

#endif

