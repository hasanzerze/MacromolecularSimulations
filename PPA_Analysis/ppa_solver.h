#ifndef PPA_SOLVER_H
#define PPA_SOLVER_H

#include <tuple>
#include <vector>
#include "vector3d.h"
#include "config.h"  // For Bond, LJParams, PPASolverParams
#include "cell_list.h" 

void computeBondForces(const std::vector<Vector3D>& positions,
                       const std::vector<config::Bond>& bonds,
                       std::vector<Vector3D>& forces,
                       const Vector3D& box_size);  // Add this

void computeLJForces_celllist(const std::vector<Vector3D>& positions,
                               std::vector<Vector3D>& forces,
                               const config::LJParams& lj,
                               const Vector3D& box_size,
                               const std::vector<int>& atom_to_chain);

void runFIREMinimizationStep(std::vector<Vector3D>& positions,
                             std::vector<Vector3D>& velocities,
                             const std::vector<config::Bond>& bonds,
                             const std::vector<std::vector<int>>& chains,
                             const config::LJParams& lj,
                             const config::PPASolverParams& params,
                             const Vector3D& box_size);

std::tuple<double, double> runPPA(
    std::vector<Vector3D>& positions,
    const std::vector<config::Bond>& bonds,
    const std::vector<std::vector<int>>& chains,
    const config::LJParams& lj,
    const config::PPASolverParams& params,
    const Vector3D& box_size,
    int frame
);

#endif

