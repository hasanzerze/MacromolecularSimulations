#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
using namespace std;

// Constants
const double LJ_EPSILON = 1.0;
const double LJ_SIGMA = 1.0;
const double sigma3 = LJ_SIGMA * LJ_SIGMA * LJ_SIGMA;
const double sigma6 = sigma3 * sigma3;
const double sigma12 = sigma6 * sigma6;
const double temperature = 0.85;
const double k_B = 1.0; // Boltzmann's constant

const double BOX_LENGTH = 8.22; // Simulation box size
const double BOX_VOLUME = BOX_LENGTH * BOX_LENGTH * BOX_LENGTH;
const int NUM_PARTICLES = 500;  // Number of particles
const double density = NUM_PARTICLES / BOX_VOLUME;
const double RCUT = 3.0 * LJ_SIGMA; // LJ cutoff
#define NUM_CELLS (static_cast<int>(BOX_LENGTH / RCUT))
const double STEP_SIZE = 0.1;  // Define the step size

struct Vector3d {
    double x, y, z;

    // Constructor
    Vector3d(double x_=0, double y_=0, double z_=0) : x(x_), y(y_), z(z_) {}

    // Addition
    Vector3d operator+(const Vector3d& v) const {
        return Vector3d(x + v.x, y + v.y, z + v.z);
    }

    // Subtraction
    Vector3d operator-(const Vector3d& v) const {
        return Vector3d(x - v.x, y - v.y, z - v.z);
    }

    // Scalar Multiplication
    Vector3d operator*(double scalar) const {
        return Vector3d(x * scalar, y * scalar, z * scalar);
    }

    // Generate a random Vector3d with values between -1 and 1
    static Vector3d Random() {
        return Vector3d(
            2.0 * (rand() / (double)RAND_MAX) - 1.0, // Random value in [-1,1]
            2.0 * (rand() / (double)RAND_MAX) - 1.0,
            2.0 * (rand() / (double)RAND_MAX) - 1.0
        );
    }

    // Generate a random Vector3d with values between -0.5 and 0.5
    static Vector3d RandomHalf() {
        return Vector3d(
            (rand() / (double)RAND_MAX) - 0.5, // Random value in [-0.5,0.5]
            (rand() / (double)RAND_MAX) - 0.5,
            (rand() / (double)RAND_MAX) - 0.5
        );
    }

    // Element-wise absolute value (cwiseAbs)
    Vector3d cwiseAbs() const {
        return Vector3d(fabs(x), fabs(y), fabs(z));
    }

};

// Particle structure
struct Particle {
    Vector3d pos; // Position as a 3D vector

    // Default constructor
    Particle() : pos(0.0, 0.0, 0.0) {}
    
    Particle(double x, double y, double z) : pos{x, y, z} {} // Constructor
//    double x, y, z;
    int cellIndex;  // Store the particle's current cell
};

// Simulation class
class Simulation {
private:
    vector<Particle> particles;
    double temperature;
    int numAcceptedMoves = 0; // Track number of accepted moves
    int numAttemptedMoves = 0; // Track number of attempted moves
    
public:
    Simulation(double temp) : temperature(temp) { initializeParticles(); }
    void initializeParticles(); // Setup initial particle positions
    vector<int> findNeighbors(int index); // Find neighbors within cutoff
    double computeLocalEnergy(int index); // Compute LJ energy for a particle
    double computeTotalEnergy();
    void performTranslationMove(); // Attempt a particle displacement
    bool metropolisCriteria(double deltaE); // Metropolis acceptance
    double computePressure(); // Compute system pressure
    void run(int numSteps); // Main simulation loop
    void printParticles(); // Function to print particle positions
    double computeVirialPressure();
    double computeLongRangeCorrectionEnergy();
    double computeLongRangeCorrectionPressure();
    void applyPeriodicBoundary(double& x, double& y, double& z);
    void applyMinimumImage(double& x, double& y, double& z);
    bool acceptMove(int index, double newX, double newY, double newZ);

    void attemptMove(int index);  // Declare the function

    vector<vector<int>> cellList; // Stores particle indices for each cell
    void buildCellList();
};

void Simulation::applyPeriodicBoundary(double& x, double& y, double& z) {
    x -= BOX_LENGTH * floor(x / BOX_LENGTH);
    y -= BOX_LENGTH * floor(y / BOX_LENGTH);
    z -= BOX_LENGTH * floor(z / BOX_LENGTH);
}

void Simulation::applyMinimumImage(double& x, double& y, double& z) {
    x -= BOX_LENGTH * round(x / BOX_LENGTH);
    y -= BOX_LENGTH * round(y / BOX_LENGTH);
    z -= BOX_LENGTH * round(z / BOX_LENGTH);
}

std::tuple<int, int, int> computeCellIndex(const Vector3d& pos, double L, int Nc) {
    double cellSize = L / Nc;

    int cx = static_cast<int>(pos.x / cellSize) % Nc;
    int cy = static_cast<int>(pos.y / cellSize) % Nc;
    int cz = static_cast<int>(pos.z / cellSize) % Nc;

    return {cx, cy, cz};
}

// Function Implementations (To be filled in later)
void Simulation::initializeParticles() {
    particles.resize(NUM_PARTICLES);
    double minSeparation = 0.8 * LJ_SIGMA; // Minimum allowed distance

    for (int i = 0; i < NUM_PARTICLES; ++i) {
        bool positionAccepted = false;
        int maxAttempts = 1000; // Prevent infinite loops
        int attempts = 0;

        while (!positionAccepted && attempts < maxAttempts) {
            // Generate a random position inside the box
            Vector3d randomVector = Vector3d::Random().cwiseAbs() * BOX_LENGTH;

            // Check for overlap
            positionAccepted = true;
            for (int j = 0; j < i; ++j) {
		
                Vector3d dr = randomVector - particles[j].pos;

                // Apply periodic boundary conditions
                applyMinimumImage(dr.x, dr.y, dr.z);

                double distanceSquared = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                if (distanceSquared < minSeparation * minSeparation) {
                    positionAccepted = false;
                    break; // Reject and retry
                }
            }

            if (positionAccepted) {
                particles[i].pos = randomVector;
            }

            attempts++;
        }

        if (attempts >= maxAttempts) {
            cerr << "Warning: Could not place particle " << i 
	    << " after " << maxAttempts << " attempts. Try reducing density.\n";
            exit(1); // Exit if placement is impossible
        }
    }

    cout << "Initialized " << NUM_PARTICLES << " particles randomly without overlaps.\n";
}

void Simulation::buildCellList() {
    // Resize cell list for 3D grid
    cellList.clear();
    cellList.resize(NUM_CELLS * NUM_CELLS * NUM_CELLS);

    // Assign particles to their respective cells
    for (size_t i = 0; i < particles.size(); ++i) {

	auto [cx, cy, cz] = computeCellIndex(particles[i].pos, BOX_LENGTH, NUM_CELLS);

        int cellIndex = cx + NUM_CELLS * (cy + NUM_CELLS * cz);
        particles[i].cellIndex = cellIndex;  // Store the computed cell index
        cellList[cellIndex].push_back(i);
    }
}

vector<int> Simulation::findNeighbors(int index) {
    vector<int> neighbors;
    double rcutSquared = RCUT * RCUT;

    // Compute cell indices for target particle
    auto [cx, cy, cz] = computeCellIndex(particles[index].pos, BOX_LENGTH, NUM_CELLS);

    // Check all neighboring cells (including current cell)
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                int nx = (cx + dx + NUM_CELLS) % NUM_CELLS;
                int ny = (cy + dy + NUM_CELLS) % NUM_CELLS;
                int nz = (cz + dz + NUM_CELLS) % NUM_CELLS;

                int neighborCell = nx + NUM_CELLS * (ny + NUM_CELLS * nz);

                // Check all particles in the neighbor cell
                for (int j : cellList[neighborCell]) {
                    if (j == index) continue; // Skip self

		    Vector3d dr = particles[j].pos - particles[index].pos;

		    applyMinimumImage(dr.x, dr.y, dr.z);

                    double distanceSquared = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                    if (distanceSquared < rcutSquared) {
                        neighbors.push_back(j);
                    }
                }
            }
        }
    }

    return neighbors;
}

void Simulation::attemptMove(int index) {

    Vector3d randomVector = (Vector3d::RandomHalf()) * STEP_SIZE;

    Particle& p = particles[index];

    // Compute new position
    Vector3d newCoords = p.pos + randomVector;

    // Apply periodic boundary conditions
    applyPeriodicBoundary(newCoords.x, newCoords.y, newCoords.z);

    // Compute new cell index
    auto [cx, cy, cz] = computeCellIndex(newCoords, BOX_LENGTH, NUM_CELLS);

    int newCellIndex = cx + NUM_CELLS * (cy + NUM_CELLS * cz);

    // Ensure the new cell index is within bounds
    if (newCellIndex < 0 || newCellIndex >= static_cast<int>(cellList.size())) {
        std::cerr << "Error: newCellIndex is out of bounds!" << std::endl;
        return;
    }

    // Accept move with Metropolis criteria
    if (acceptMove(index, newCoords.x, newCoords.y, newCoords.z)) {
        // Remove particle from old cell
        auto& oldCell = cellList[p.cellIndex];
        oldCell.erase(std::remove(oldCell.begin(), oldCell.end(), index), oldCell.end());

        // Update particle's position
        p.pos = newCoords;
        p.cellIndex = newCellIndex;

        // Add particle to new cell
        cellList[newCellIndex].push_back(index);
    }
}

bool Simulation::acceptMove(int index, double newX, double newY, double newZ) {
    double oldEnergy = computeLocalEnergy(index); // Energy at current position

    // Temporarily move the particle
    Vector3d oldPos = particles[index].pos;
    particles[index].pos = {newX, newY, newZ};

    double newEnergy = computeLocalEnergy(index); // New energy after move

    // Revert position (we only move it if accepted)
    particles[index].pos = oldPos;

    double deltaE = newEnergy - oldEnergy;

    // Increment attempted move count
    numAttemptedMoves++;

    if (deltaE < 0) {
        numAcceptedMoves++; // Increment accepted move count
        return true; // Always accept if energy decreases
    } else {
        double p = exp(-deltaE / temperature);
        if ((rand() / (double)RAND_MAX) < p) {
            numAcceptedMoves++; // Increment accepted move count
            return true; // Accept with probability exp(-deltaE/kT)
        }
    }
    return false; // Reject the move
}

double Simulation::computeLocalEnergy(int index) {
    const Particle& p = particles[index];
    double energy = 0.0;

    // Get relevant neighbors from cell lists (optimized)
    vector<int> neighbors = findNeighbors(index);

    for (int neighborIndex : neighbors) {
            if (neighborIndex != index) {
                const Particle& neighbor = particles[neighborIndex];

                // Compute the distance between particles p and neighbor
                Vector3d dr = p.pos - neighbor.pos;

                // Apply minimum image conditions
                applyMinimumImage(dr.x, dr.y, dr.z);

                // Compute the squared distance
                double r2 = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;

                // If within cutoff radius
                if (r2 < RCUT * RCUT) {
                    double r6 = r2 * r2 * r2;
                    double r12 = r6 * r6;

                    // Lennard-Jones potential
                    energy += 4. * LJ_EPSILON * (sigma12/r12 - sigma6/r6);
                }
            }
        }
    return energy;
}

double Simulation::computeTotalEnergy() {
    double totalEnergy = 0.0;

    // Iterate over all particles
    for (int i = 0; i < NUM_PARTICLES; ++i) {
    	totalEnergy += computeLocalEnergy(i);
    }

    // Apply long-range correction
    totalEnergy += computeLongRangeCorrectionEnergy();

    return totalEnergy / 2.;	// We divided by 2 to avoid double counting of pairs
    }

void Simulation::run(int numSteps) {
    cout << "Running Monte Carlo simulation for " << numSteps << " steps...\n";
    double Psum = 0.;
    int nsample = 0;
    // Main Metropolis loop
    for (int step = 0; step < numSteps; ++step) {
        // Choose a random particle to move
        int index = rand() % NUM_PARTICLES;

        // Attempt to move the particle
        attemptMove(index);
        // Optionally, you can compute and print thermodynamic properties:
        if (step % 10000 == 0) {
		double totalEnergy = computeTotalEnergy();
		double pressure = computeVirialPressure();
        	double acceptanceRatio = static_cast<double>(numAcceptedMoves) / numAttemptedMoves;
		Psum += pressure;
		nsample += 1;

        	cout << "Step " << step << " - Total Energy: " << totalEnergy 
		     << ", Acceptance Ratio = " << acceptanceRatio
		     << ", Pressure = " << pressure << endl;
        }
    }

    cout << "Monte Carlo simulation finished." << endl;
    cout << "Average pressure = " << Psum / nsample << endl;
}

double Simulation::computeLongRangeCorrectionEnergy() {
    double rc3 = pow(RCUT, 3);
    double rc9 = rc3 * rc3 * rc3;

    double correction = (8.0 * M_PI * NUM_PARTICLES * density * LJ_EPSILON * pow(LJ_SIGMA, 3) / 3.0) *
                        ((1.0 / 3.0) * (pow(LJ_SIGMA, 9) / rc9) - (pow(LJ_SIGMA, 3) / rc3));

    return correction;
}

double Simulation::computeVirialPressure() {
    double virialSum = 0.0;
    
    for (int i = 0; i < NUM_PARTICLES; ++i) {
        for (int j = i + 1; j < NUM_PARTICLES; ++j) {  // Avoid double counting
	    Vector3d dr = particles[i].pos - particles[j].pos;

            // Apply periodic boundary conditions (assuming cubic box)
            applyMinimumImage(dr.x, dr.y, dr.z);

            double r2 = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
            if (r2 < RCUT * RCUT) {  // Apply LJ cutoff
                double r2_inv = 1.0 / r2;
                double r6_inv = sigma6 * r2_inv * r2_inv * r2_inv;
                double force_mag = 24.0 * LJ_EPSILON * (2.0 * r6_inv * r6_inv - r6_inv) * r2_inv;
                
                // Compute force vector
                double fx = force_mag * dr.x;
                double fy = force_mag * dr.y;
                double fz = force_mag * dr.z;

                // Accumulate virial sum: Fij.rij
                virialSum += fx * dr.x + fy * dr.y + fz * dr.z;
            }
        }
    }
    
    // Compute pressure using the virial equation
    double idealGasTerm = (NUM_PARTICLES * k_B * temperature) / (BOX_VOLUME);
    double interactionTerm = virialSum / (3.0 * BOX_VOLUME);
    
    return idealGasTerm + interactionTerm + computeLongRangeCorrectionPressure();
}

double Simulation::computeLongRangeCorrectionPressure() {
    double rc3 = pow(RCUT, 3);
    double rc9 = rc3 * rc3 * rc3;

    double correction = (16.0 * M_PI * density * density * LJ_EPSILON * pow(LJ_SIGMA, 3) / 3.0) *
                        ((2.0 / 3.0) * (pow(LJ_SIGMA, 9) / rc9) - (pow(LJ_SIGMA, 3) / rc3));

    return correction;
}

void Simulation::printParticles() {
    cout << "Particle Positions:\n";
    for (size_t i = 0; i < particles.size(); ++i) {
        cout << "Particle " << i << ": " << particles[i].pos.x << ", "
	     << particles[i].pos.y << ", " << particles[i].pos.z << "\n";
    }
}

// Main function
int main() {
    srand(time(0)); // Seed random number generator

    Simulation mcSim(temperature);

    mcSim.printParticles(); // Print initial positions

    mcSim.buildCellList(); // Build cell list before neighbor search

    mcSim.run(1000000); // Run simulation for 1,000,000 steps
    
    return 0;
}


