#include "ppa_solver.h"
#include <iostream>
#include <cmath>
#include "cell_list.h"
#include <tuple>

extern "C" {
    #include <xdrfile_xtc.h>
}

// Compute harmonic bond force between bonded monomers
void computeBondForces(const std::vector<Vector3D>& positions,
                       const std::vector<config::Bond>& bonds, double k_bond,
                       std::vector<Vector3D>& forces, const Vector3D& box_size) {
    for (const auto& bond : bonds) {
        Vector3D r = positions[bond.j] - positions[bond.i];
	r = apply_minimum_image(r, box_size);
        double dist = r.norm();
	Vector3D f = r * ( - k_bond );
        forces[bond.i] -= f;
        forces[bond.j] += f;
    }
}

void computeLJForces_celllist(const std::vector<Vector3D>& positions,
                               std::vector<Vector3D>& forces,
                               const config::LJParams& lj,
                               const Vector3D& box_size,
			       const std::vector<int>& atom_to_chain) {
    const double sigma2 = lj.sigma * lj.sigma;

    double cutoff = std::pow(2.0, 1.0 / 6.0) * lj.sigma;
    double cutoff_sq = cutoff * cutoff;

    CellList grid(cutoff, box_size);
    grid.build(positions, box_size);

    for (const auto& [cellIdx, atomList] : grid.cells) {
        auto neighbors = grid.getNeighboringCells(cellIdx);
        for (int i : atomList) {
            for (const CellIndex& nbr : neighbors) {
		const auto& neighbor_atoms = grid.cells[nbr]; // will default to empty if missing
		for (int j : neighbor_atoms) {
                    if (j <= i) continue;  // avoid double-count

		    // Skip if same chain
                    if (atom_to_chain[i] == atom_to_chain[j]) continue;

                    Vector3D r = positions[i] - positions[j];
                    r = apply_minimum_image(r, box_size);
                    double r2 = r.dot(r);

                    if (r2 < cutoff_sq) {
                        double inv_r2 = sigma2 / r2;
                        double inv_r6 = inv_r2 * inv_r2 * inv_r2;
                        double fmag = 48.0 * lj.epsilon * inv_r6 * (inv_r6 - 0.5) / r2;
                        Vector3D f = r * fmag;
                        forces[i] += f;
                        forces[j] -= f;
                    }
                }
            }
        }
    }
}

void runFIREMinimizationStep(std::vector<Vector3D>& positions,
                             std::vector<Vector3D>& velocities,
                             const std::vector<config::Bond>& bonds,
                             const std::vector<std::vector<int>>& chains,
                             const config::LJParams& lj,
                             const config::PPASolverParams& params,
                             const Vector3D& box_size) {
    size_t N = positions.size();
    size_t Nchain = chains.size();
    size_t chainsize = chains[0].size();
    std::vector<Vector3D> forces(N, Vector3D{0.0, 0.0, 0.0});
    std::fill(velocities.begin(), velocities.end(), Vector3D{0.0, 0.0, 0.0});

    double dt = params.timestep;
    int maxiter = params.max_iterations;
    const double tol = params.bond_length_tolerance;
    const double force_tol = 1e-3;
    const double force_tol2 = force_tol * force_tol;

    double alpha = 0.1;
    double dt_max = dt;
    const double finc = 1.1, fdec = 0.5, falpha = 0.99;
    int n_since_negative = 0;
    int n_min = 5;

    std::vector<int> atom_to_chain(N, -1);
    for (size_t c = 0; c < Nchain; ++c) {
        for (int atom : chains[c]) {
             atom_to_chain[atom] = static_cast<int>(c);
        }
    }

    for (int iter = 0; iter < maxiter; ++iter) {
        std::fill(forces.begin(), forces.end(), Vector3D{0, 0, 0});
        computeBondForces(positions, bonds, params.ppa_bond_const, forces, box_size);
        computeLJForces_celllist(positions, forces, lj, box_size, atom_to_chain);

        double max_force_norm = 0.0;
        double P = 0.0;

        for (size_t c = 0; c < Nchain; ++c) {
            const std::vector<int>& chain = chains[c];
            for (size_t m = 0; m < chainsize; ++m) {
                int i = chain[m];
                // Skip if it's the first or last monomer in the chain
                if (m == 0 || m == chainsize - 1) continue;
                velocities[i] += forces[i] * dt;
                positions[i] += velocities[i] * dt;

                P += velocities[i].dot(forces[i]);
                double f2 = forces[i].dot(forces[i]);
                if (f2 > max_force_norm) max_force_norm = f2;
            }
        }
//        std::cout << "[DEBUG] Iter " << iter
//                  << " | maxF = " << std::sqrt(max_force_norm)
//                  << " | dt = " << dt
//                  << " | alpha = " << alpha
//                  << " | P = " << P << "\n";

        if (P > 0.0) {
                n_since_negative++;
                if (n_since_negative > n_min) {
                    dt = std::min(dt * finc, dt_max);
                    alpha *= falpha;
                }
        } else {
                dt *= fdec;
                alpha = 0.1;
                n_since_negative = 0;

                // Zero velocities
                for (size_t i = 0; i < velocities.size(); ++i) {
                    velocities[i] = Vector3D{0, 0, 0};
                }
        }

        for (size_t c = 0; c < Nchain; ++c) {
            const auto& chain = chains[c];
            for (size_t m = 1; m + 1 < chainsize; ++m) {
                int i = chain[m];
                velocities[i] = velocities[i] * (1 - alpha) + forces[i].normalized() * alpha * velocities[i].norm();
            }
        }

        if (max_force_norm < force_tol2) {
            std::cout << "Converged after " << iter << " FIRE iterations (force criterion).\n";
            break;
        }
    }
}

std::tuple<double, double> runPPA(std::vector<Vector3D>& positions,
              const std::vector<config::Bond>& bonds,
              const std::vector<std::vector<int>>& chains,
              const config::LJParams& lj,
              const config::PPASolverParams& params,
	      const Vector3D& box_size, int frame) {

    const int max_rounds = 100;
    const double Lpp_tol = 1e-3;

    double previous_Lpp = -1.0;
    double current_Lpp = -1.0;
    double Ree2_av = -1.0;

    // Open XTC file to write multiple frames
    std::string xtcfilename = "ppa_" + std::to_string(frame) + ".xtc";
    XDRFILE* xtc_out = xdrfile_open(xtcfilename.c_str(), "w");
    if (!xtc_out) {
        std::cerr << "Could not open XTC file for writing.\n";
        return {-1.0, -1.0};
    }

    matrix box_matrix;
    box_matrix[0][0] = box_size.x; box_matrix[0][1] = 0; box_matrix[0][2] = 0;
    box_matrix[1][0] = 0; box_matrix[1][1] = box_size.y; box_matrix[1][2] = 0;
    box_matrix[2][0] = 0; box_matrix[2][1] = 0; box_matrix[2][2] = box_size.z;

    rvec* x = new rvec[positions.size()];

    for (int round = 0; round < max_rounds; ++round) {
        std::cout << "\n=== PPA ROUND " << round + 1 << " ===\n";

        // Zero velocities before each round (if FIRE is used)
	std::vector<Vector3D> velocities(positions.size(), Vector3D{0.0, 0.0, 0.0});

	// MINIMIZATION LOOP (FIRE, SD, etc.)
	runFIREMinimizationStep(positions, velocities, bonds, chains, lj, params, box_size);

    // Compute average contour length per chain (entanglement length proxy)
	double totalContourLength = 0.0;
	int totalBonds = 0;
	double totalendtoend2 = 0.0;

	for (const auto& chain : chains) {
		Vector3D re = positions[chain.back()] - positions[chain.front()];
		re = apply_minimum_image(re, box_size);
		totalendtoend2 += re.dot(re);
        	for (size_t i = 1; i < chain.size(); ++i) {
        	    Vector3D diff = positions[chain[i]] - positions[chain[i - 1]];
		    diff = apply_minimum_image(diff, box_size);
        	    totalContourLength += diff.norm();
        	    totalBonds++;
        	}
    	}

	Ree2_av = totalendtoend2 / chains.size();
    	current_Lpp = (totalBonds > 0) ? totalContourLength / chains.size() : -1.0;

    	std::cout << "  --> Estimated Lpp: " << current_Lpp << " Mean square end-to-end distance: " << Ree2_av << "\n";

        // Convergence Check
        if (round > 0 && std::abs(current_Lpp - previous_Lpp)/current_Lpp < Lpp_tol) {
            std::cout << " Lpp converged after " << round + 1 << " rounds.\n";
            break;
        }
        previous_Lpp = current_Lpp;

	// Fill rvec for writing
    	for (size_t i = 0; i < positions.size(); ++i) {
        	x[i][0] = positions[i].x;
        	x[i][1] = positions[i].y;
        	x[i][2] = positions[i].z;
    	}

    	// Write one frame per round
	write_xtc(xtc_out, positions.size(), round, 0.0f, box_matrix, x, 1000.0f);
    }

    delete[] x;
    xdrfile_close(xtc_out);

    return {current_Lpp, Ree2_av};
}


