#include <cmath>
#include <iostream>
#include <xdrfile.h>
#include <xdrfile_xtc.h>

#include <complex>
#include <vector>
#include <stdexcept>
#include <algorithm>

using namespace std;

#define PI 3.14159265

// Function to compute double factorial (n!!)
int doubleFactorial(int n) {
    if (n <= 0) return 1;
    int result = 1;
    for (int i = n; i > 1; i -= 2) {
        result *= i;
    }
    return result;
}

double associatedLegendrePolynomial(int l, int m, double x) {
    if (m < 0 || m > l || std::abs(x) > 1.0) {
        throw std::invalid_argument("Invalid arguments for associated Legendre polynomial.");
    }

    // Create a 2D vector to store P_l^m(x)
    std::vector<std::vector<double>> P(l + 1, std::vector<double>(l + 1, 0.0));

    // Base cases
    P[0][0] = 1.0;
    if (l > 0) {
        P[1][0] = x;
    }
    for (int k = 0; k <= l; ++k) {
        P[k][k] = std::pow(-1.0, k) * doubleFactorial(2 * k - 1) * std::pow(1 - x * x, k / 2.0);
    }

    // Fill the table using the recurrence relation
    for (int n = 2; n <= l; ++n) {
        for (int k = 0; k < n; ++k) {
            P[n][k] = ((2 * n - 1) * x * P[n - 1][k] - (n + k - 1) * P[n - 2][k]) / (n - k);
        }
    }

    return P[l][m];
}

// Function to compute the spherical harmonics Y_l^m(theta, phi)
std::complex<double> sphericalHarmonic(int l, int m, double theta, double phi) {
    // Compute the associated Legendre polynomial P_l^m(cos(theta))
    //double P_lm = associatedLegendrePolynomial(l, std::abs(m), std::cos(theta));
    double P_lm = assoc_legendre(l, std::abs(m), std::cos(theta));

    // Compute the normalization factor
    float sign = 1.0;
    if ( m != 0 ) {
        sign = pow((-abs(m)/m), abs(m) );
    }

    double normalization = sign * std::sqrt((2.0 * l + 1.0) / (4.0 * M_PI) * std::tgamma(l - std::abs(m) + 1) / std::tgamma(l + std::abs(m) + 1));

    // Compute the complex exponential factor
    std::complex<double> expFactor = std::exp(std::complex<double>(0, m * phi));

    // Combine everything to get Y_l^m(theta, phi)
    return normalization * P_lm * expFactor;
}

// Function to compute a vector of spherical harmonics Y_l^m for m = -l to l
std::vector<std::complex<double>> sphericalHarmonicsVector(int l, double theta, double phi) {
    std::vector<std::complex<double>> Y_lm_vector;
    for (int m = -l; m <= l; ++m) {
        Y_lm_vector.push_back(sphericalHarmonic(l, m, theta, phi));
    }
    return Y_lm_vector;
}

// Function to sum two vectors element-wise
std::vector<std::complex<double>> sumVectors(const std::vector<std::complex<double>>& vec1, const std::vector<std::complex<double>>& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Vectors must be of the same size to sum element-wise.");
    }
    
    std::vector<std::complex<double>> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <xtc file>" << std::endl;
        return 1;
    }

    const char *filename = argv[1];
    XDRFILE *xd;
    int natoms;
    int step;
    float time;
    float xij;
    float yij;
    float zij;
    float rxy;
    float rij;
    float dist2;
    float r_neigh = 1.4;
    float r_neigh2 = r_neigh * r_neigh;
    float theta;
    float phi;
    int l = 4;
    //std::complex<double> Y_lm;
    std::vector<std::complex<double>> Y_lm_vector(2 * l + 1, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> new_Y_lm_vector(2 * l + 1, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> q_lm(2 * l + 1, std::complex<double>(0.0, 0.0));

    matrix box;
    rvec *x = NULL;
    float prec;

    // Open the XTC file
    xd = xdrfile_open(filename, "r");
    if (xd == NULL) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return 1;
    }

    // Read the number of atoms
    if (read_xtc_natoms(const_cast<char*>(filename), &natoms) != exdrOK) {
    	std::cerr << "Error: Unable to read number of atoms from " << filename << std::endl;
    	xdrfile_close(xd);
    	return 1;
    }

	// Allocate memory for coordinates
    x = (rvec *)malloc(natoms * sizeof(rvec));
    if (x == NULL) {
        std::cerr << "Error: Unable to allocate memory for coordinates" << std::endl;
        xdrfile_close(xd);
        return 1;
    }

    // Read frames
    while (read_xtc(xd, natoms, &step, &time, box, x, &prec) == exdrOK) {
        std::cout << "Time: " << time << " ps" << " Box size = " << box[0][0] << ", "<< box[1][1] << ", " << box[2][2] << ", " << std::endl;
        for (int i = 0; i < natoms; ++i) {
            //std::cout << "Atom " << i + 1 << ": (" << x[i][0] << ", " << x[i][1] << ", " << x[i][2] << ")" << std::endl;
            std::fill(Y_lm_vector.begin(), Y_lm_vector.end(), std::complex<double>(0.0, 0.0));
	    int n_neigh = 0;
	    for (int j = 0; j < natoms; ++j) {
		xij = x[j][0] - x[i][0];
		yij = x[j][1] - x[i][1];
		zij = x[j][2] - x[i][2];
		xij = xij - box[0][0] * round(xij / box[0][0]);
		yij = yij - box[1][1] * round(yij / box[1][1]);
		zij = zij - box[2][2] * round(zij / box[2][2]);

		dist2 = xij * xij + yij * yij + zij * zij;
		if (i != j and r_neigh2 > dist2 ) {
			n_neigh += 1;
			rij = sqrt(dist2);
			rxy = sqrt(xij * xij + yij * yij);
			theta = PI / 2.0;
			if ( zij != 0.0 ) {
				if ( zij > 0.0 ) {
					theta = atan(rxy / zij);
					//cout << "theta = " << theta;
				}
				if ( zij < 0.0 ) {
					theta = PI - atan( - rxy / zij);
					//cout << "theta = " << theta;
				}
			}
			if ( xij >= 0 ) {
				if ( yij >= 0 ) {
					phi = PI / 2.0;
				}
				else {
					phi = - PI / 2.0;
				}

				if ( xij != 0.0 ) {
					phi = atan ( yij / xij );
					//cout << "phi = " << phi;
				}
			}
			if ( xij < 0 ) {
				phi = PI + atan ( yij / xij );
				//cout << "phi = " << phi;
			}
			for (int m = -l; m <= l; ++m) {
				new_Y_lm_vector[m+l] = sphericalHarmonicsVector(l, theta, phi)[m+l];
			}
			// Sum the vectors element-wise
			Y_lm_vector = sumVectors(Y_lm_vector, new_Y_lm_vector);
		}
	    }
	    for (int m = -l; m <= l; ++m) {
	    	q_lm[m+l] = Y_lm_vector[m+l] / static_cast<double>(n_neigh);
		//cout << "q_4," << m << " = " << q_lm[m+l] << endl;
	    }
	    float q_l = 0;
	    for (int m = -l; m <= l; ++m) {
		q_l += norm ( q_lm[ m + l ] );
	    }
	    q_l = sqrt( 4 * PI / (2 * l + 1) * q_l);
	    cout << "ql = " << q_l << endl;
	    //cout << " Number of neighbors = " << n_neigh << endl;
        }
    }

    // Clean up
    free(x);
    xdrfile_close(xd);

    return 0;
}

