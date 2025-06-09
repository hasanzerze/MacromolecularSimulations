#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;

// Global parameters (initialized later)
int N;  // Number of spins
double J;  // Interaction strength
double h;  // External field
double T;  // Temperature
int MC_STEPS;  // Monte Carlo steps
double invT;         // 1 / Temperature

int accepted_flips = 0;  // Counter for accepted flips
double Esum = 0;
double E2sum = 0;
int datafreq = 10;
int ndata = 0;

// Random number generator
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<double> uniform(0.0, 1.0);
vector<int> spins;

// Function to read parameters from a file
bool read_parameters(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error: Could not open " << filename << endl;
        return false;
    }

    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string param;
        char eq;
        if (ss >> param >> eq) {
            if (eq != '=')
                continue;

            if (param == "N") ss >> N;
            else if (param == "J") ss >> J;
            else if (param == "h") ss >> h;
            else if (param == "T") ss >> T;
            else if (param == "MC_STEPS") ss >> MC_STEPS;
        }
    }
    file.close();

    invT = 1.0 / T; // Compute inverse temperature
    return true;
}

// Initialize spins randomly
void initialize_spins() {
    spins.resize(N);
    for (int i = 0; i < N; i++) {
        spins[i] = (uniform(gen) < 0.5) ? -1 : 1;
    }
}

// Compute total energy of the system
double compute_energy() {
    double E = 0.0;
    for (int i = 0; i < N; i++) {
        int next = (i + 1) % N;  // Periodic boundary condition
        E -= J * spins[i] * spins[next];
    }
    return E - h * accumulate(spins.begin(), spins.end(), 0.0);
}

// Perform a Metropolis step
void metropolis_step() {
    uniform_int_distribution<int> rand_site(0, N - 1);

    for (int step = 0; step < N; step++) {  // N trials per MC step
        int i = rand_site(gen);  // Choose a random spin
        int left = (i - 1 + N) % N;
        int right = (i + 1) % N;
        double dE = 2 * J * spins[i] * (spins[left] + spins[right]) + 2 * h * spins[i];

        if (dE <= 0 || uniform(gen) < exp(-invT * dE)) {
            spins[i] *= -1;  // Flip the spin
	    accepted_flips++; // Count the accepted move
        }
    }
}

// Run the Monte Carlo simulation
void run_simulation() {
    initialize_spins();

    // Thermalization
    for (int i = 0; i < MC_STEPS / 2; i++) {
        metropolis_step();
    }

    // Measurement
    ofstream file("ising_output.dat");
    for (int i = 0; i < MC_STEPS; i++) {
        metropolis_step();
        double E = compute_energy();
        int M = accumulate(spins.begin(), spins.end(), 0);

        file << i << " " << E << " " << M << endl;  // Save data
	if ( (i % datafreq) == 0 ) {
		Esum += E;
		E2sum += E * E;
		ndata += 1;
	}
    }
    file.close();
}

int main() {
    if (!read_parameters("input.txt")) {
        return 1;  // Exit if input file not found
    }
    run_simulation();
    cout << "Simulation complete! Data saved to ising_output.dat" << endl;
    cout << "Acceptance fraction = " << 1.0 * accepted_flips / MC_STEPS / N << endl;
    cout << "Number of accepted attempts = " << accepted_flips << endl;
    cout << "Number of MC steps = " << MC_STEPS << endl;
    cout << "<E> / N = " << Esum / ndata / N << endl;
    cout << "<E2> = " << E2sum / ndata << endl;
    cout << "Cv / N = " << invT * invT * ( ( E2sum / ndata ) - ( Esum / ndata )*( Esum / ndata ) ) / N << endl;
    return 0;
}



