#ifndef XTC_READER_H
#define XTC_READER_H

#include <vector>
#include "vector3d.h"

extern "C" {
    #include <xdrfile_xtc.h>
}

// Reads one frame and stores atom positions
bool read_xtc_frame(XDRFILE* xtc, std::vector<Vector3D>& positions, int natoms, float& time, int& step, float& prec, Vector3D& box_size);

#endif

