#include <cstring> // for strcpy
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "xtc_reader.h"
#include "vector3d.h"
#include "config.h"
#include "ppa_solver.h"

extern "C" {
    #include <xdrfile_xtc.h>
}

int main() {
    // Load simulation config (e.g., filenames, number of atoms, bonds, etc.)
    //config::load_config_from_file("input.dat");
    config::loadConfigFromFile("input.dat");

    // Open the XTC file
    XDRFILE *xtc = xdrfile_open(config::XTC_FILE.c_str(), "r");

    int natoms_xtc = 0;
    int rc = read_xtc_natoms(const_cast<char*>(config::XTC_FILE.c_str()), &natoms_xtc);

    if (rc != exdrOK) {
        std::cerr << "[ERROR] Could not read natoms from XTC file.\n";
        return 1;
    }

    std::cout << "[INFO] Atoms in XTC file: " << natoms_xtc << "\n";
    if (natoms_xtc != config::NUM_ATOMS) {
        std::cerr << "[ERROR] Mismatch: config::NUM_ATOMS = " << config::NUM_ATOMS
                  << ", but XTC contains " << natoms_xtc << "\n";
        return 1;
    }

    if (!xtc) {
        std::cerr << "Error: Cannot open XTC file.\n";
        return 1;
    }

    // Prepare storage for Lpp values
    std::vector<double> Lpp_values;
    std::vector<double> Ne_values;
    std::vector<double> app_values;
    std::vector<double> bpp_values;
    std::vector<double> Ree2_values;

    int natoms = config::NUM_ATOMS;
    int natoms_per_chain = config::MONOMERS_PER_CHAIN;

    std::vector<Vector3D> positions(natoms);
    float time;
    int step;
    float prec = 0.001;

    int frameIndex = 0;

    std::ofstream fout("Lpp_vs_time.csv");
    fout << "#Frame, Time, Lpp, app, bpp, Ne\n";

    Vector3D box_size;

    while (read_xtc_frame(xtc, positions, natoms, time, step, prec, box_size)) {
        if (frameIndex % config::FRAME_STRIDE != 0) {
            ++frameIndex;
            continue;
        }

        std::cout << "Processing frame " << frameIndex << " at time " << time << " ps...\n";

        auto [Lpp, Ree2] = runPPA(positions, config::BONDS, config::CHAINS, config::LJ_PARAMS, config::PPA_PARAMS, box_size, frameIndex);
	double app = Ree2 / Lpp;
	double bpp = Lpp /  natoms_per_chain;
	double Ne = app / bpp;
        if (Lpp > 0.0) {
            Lpp_values.push_back(Lpp);
	    Ree2_values.push_back(Ree2);
            app_values.push_back(app);
            bpp_values.push_back(bpp);
            Ne_values.push_back(Ne);
        }

	fout << frameIndex << "," << time << "," << Lpp << "," << app <<"," << bpp <<"," << Ne << "\n";

        ++frameIndex;
    }

    fout.close();

    xdrfile_close(xtc);

    // Compute average and standard deviation of Lpp
    if (!Lpp_values.empty()) {
        double sum = 0.0, sumsq = 0.0;
        for (double ne : Lpp_values) {
            sum += ne;
            sumsq += ne * ne;
        }
        double Lpp_avg = sum / Lpp_values.size();
        double Lpp_stddev = std::sqrt((sumsq / Lpp_values.size()) - Lpp_avg * Lpp_avg);

        sum = 0.0, sumsq = 0.0;
        for (double app : app_values) {
            sum += app;
            sumsq += app * app;
        }
        double app_avg = sum / app_values.size();
        double app_stddev = std::sqrt((sumsq / app_values.size()) - app_avg * app_avg);

        sum = 0.0, sumsq = 0.0;
        for (double bpp : bpp_values) {
            sum += bpp;
            sumsq += bpp * bpp;
        }
        double bpp_avg = sum / bpp_values.size();
        double bpp_stddev = std::sqrt((sumsq / bpp_values.size()) - bpp_avg * bpp_avg);

        sum = 0.0, sumsq = 0.0;
        for (double Ne : Ne_values) {
            sum += Ne;
            sumsq += Ne * Ne;
        }
        double Ne_avg = sum / Ne_values.size();
        double Ne_stddev = std::sqrt((sumsq / Ne_values.size()) - Ne_avg * Ne_avg);

        std::cout << "\nAverage Lpp over " << Lpp_values.size() << " frames: " << Lpp_avg
                  << " ± " << Lpp_stddev 
		  << "\nAverage app over " << app_values.size() << " frames: " << app_avg
                  << " ± " << app_stddev
                  << "\nAverage bpp over " << bpp_values.size() << " frames: " << bpp_avg
                  << " ± " << bpp_stddev
                  << "\nAverage Ne over " << Ne_values.size() << " frames: " << Ne_avg
                  << " ± " << Ne_stddev
		  << std::endl;
    } else {
        std::cout << "No valid frames were processed.\n";
    }

    return 0;
}

