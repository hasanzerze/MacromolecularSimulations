#include "xtc_reader.h"
#include <iostream>

extern "C" {
    #include <xdrfile_xtc.h>
}

bool read_xtc_frame(XDRFILE* xtc, std::vector<Vector3D>& positions, int natoms, float& time, int& step, float& prec, Vector3D& box_size) {
    if (!xtc) {
        std::cerr << "[ERROR] XTC file pointer is null.\n";
        return false;
    }

    if (positions.size() != static_cast<size_t>(natoms)) {
        std::cerr << "[WARNING] Resizing positions array to match natoms = " << natoms << "\n";
        positions.resize(natoms);
    }

    rvec* x = new rvec[natoms];
    matrix box;
    int result = read_xtc(xtc, natoms, &step, &time, box, x, &prec);

    if (result != exdrOK) {
        delete[] x;
        return false;
    }

    for (int i = 0; i < natoms; ++i) {
        positions[i] = Vector3D(x[i][0], x[i][1], x[i][2]);
    }

    box_size.x = box[0][0];
    box_size.y = box[1][1];
    box_size.z = box[2][2];

    delete[] x;
    return true;
}

