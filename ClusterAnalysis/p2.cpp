#include <iostream>
#include <map>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <xdrfile.h>
#include <xdrfile_xtc.h>
#include <stack>
#include <iomanip>
#include <algorithm>

using namespace std;

// Iterative DFS to find all atoms in a cluster
void dfs_iterative(int start,
                   const vector<vector<int>>& adj,
                   vector<bool>& visited,
                   vector<int>& cluster) {
    stack<int> s;
    s.push(start);

    while (!s.empty()) {
        int current = s.top();
        s.pop();

        if (visited[current]) continue;

        visited[current] = true;
        cluster.push_back(current);

        for (int neighbor : adj[current]) {
            if (!visited[neighbor]) {
                s.push(neighbor);
            }
        }
    }
}

void printHistogram(const std::vector<double>& data, double binSize) {
    // Create a map to hold the histogram
    std::map<int, int> histogram;

    // Count occurrences in each bin
    for (double number : data) {
         // Determine which bin the number belongs to
         int bin = static_cast<int>(std::floor(number / binSize));
         histogram[bin]++;
    }
    // Output file for writing p2 histogram:
    std::ofstream outputfile_hist;
    outputfile_hist.open("hist.out");
    outputfile_hist << "#p2     Counts" << "\n";

    // Print the histogram
    for (const auto& pair : histogram) {
         outputfile_hist << (pair.first * binSize) << " " << pair.second << "\n";
    }

    outputfile_hist.close();
}

struct Bond {
    int atom_index;   // Atom involved in the bond
    float p2;
    float mag;
    bool crystal;
    rvec bond_vector; // 3D vector associated with the bond (rvec is defined as float[3])
    std::vector<int> neighbors;  // Indices of neighboring bonds
};

// Function to calculate the Euclidean distance using the Minimum Image Convention (MIC)
float calculate_distance2_MIC(const rvec &atom1, const rvec &atom2, const matrix box) {
     float dx = atom2[0] - atom1[0];
     float dy = atom2[1] - atom1[1];
     float dz = atom2[2] - atom1[2];
     // Apply Minimum Image Convention (MIC) for each dimension
     // Wrapping distances to ensure they are within the minimum image

     dx -= round(dx / box[0][0]) * box[0][0];  // Adjust x-coordinate
     dy -= round(dy / box[1][1]) * box[1][1];  // Adjust y-coordinate
     dz -= round(dz / box[2][2]) * box[2][2];  // Adjust z-coordinate

     return (dx*dx + dy*dy + dz*dz);                                             
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input file>" << std::endl;
        return 1;
    }

    // Output file for writing time dependent crystal fraction 
    std::ofstream outputfile;
    outputfile.open("crystallinity.out");
    outputfile << "# Time (ps) CrystalFraction";

    // Output XYZ file for labeling clusters:
    std::ofstream xyz_file("clustered.xyz");

    const char *input_filename = argv[1];
    std::ifstream input_file(input_filename);
    if (!input_file) {
	std::cerr << "Error: Unable to open input file " << input_filename << std::endl;
	return 1;
    }

    std::string line;
    std::string xtc_filename;
    int imol = 0;
    int iatom = 0;
    int jatom = 0;
    int monomers_per_molecule = 0;
    int num_bonds = 0;
    int nmolecules = 0;
    float rcut = 0.;
    float clustcutoff = 0.;
    float p2thresh = 0.0;
    int frame_skip = 1;
    double binSize = 0.01; // Define the size of each bin for histogram

    // Read the input file for the needed variables:

    while (std::getline(input_file, line)) {
	if (line.find("xtc_file") != std::string::npos) {
		xtc_filename = line.substr(line.find("=")+2); // Extract value after '='
	} else if (line.find("monomers_per_molecule") != std::string::npos) {
		monomers_per_molecule = std::stoi(line.substr(line.find("=")+2));
	} else if (line.find("rcut") != std::string::npos) {
		rcut = std::stof(line.substr(line.find("=")+2));
	} else if (line.find("clustcutoff") != std::string::npos) {
		clustcutoff = std::stof(line.substr(line.find("=")+2));
	} else if (line.find("p2thresh") != std::string::npos) {
		p2thresh = std::stof(line.substr(line.find("=")+2));
	} else if ( line.find("frame_skip") != std::string::npos) {
		frame_skip = std::stoi(line.substr(line.find("=")+2));
	}
    }

    float rcut2 = rcut * rcut;
    float clustcutoff2 = clustcutoff * clustcutoff;

    input_file.close();

    // Check if the variables were successfully read
    if (xtc_filename.empty() || monomers_per_molecule == 0) {
	std::cerr << "Error: Missing xtc_file or monomers_per_molecule in the input file" << std::endl;
	return 1;
    }

    std::cout << "XTC File: " << xtc_filename << std::endl;
    std::cout << "Monomers per molecule: " << monomers_per_molecule << std::endl;

    XDRFILE *xd;
    int natoms;
    int step;
    float time;
    matrix box;
    rvec *x = NULL;
    float prec;

    // Open the XTC file
    xd = xdrfile_open(const_cast<char*>(xtc_filename.c_str()), "r");
    if (xd == NULL) {
        std::cerr << "Error: Unable to open file " << xtc_filename << std::endl;
        return 1;
    }

    // Read the number of atoms
    if (read_xtc_natoms(const_cast<char*>(xtc_filename.c_str()), &natoms) != exdrOK) {
    	std::cerr << "Error: Unable to read number of atoms from " << xtc_filename << std::endl;
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

    nmolecules = natoms / monomers_per_molecule;
    num_bonds = ( monomers_per_molecule - 2 ) * nmolecules;

    std::cout << "XTC File: " << xtc_filename << std::endl;
    std::cout << "Monomers per molecule: " << monomers_per_molecule << std::endl;
    std::cout << "Number of molecules: " << nmolecules << std::endl;
    std::cout << "Number of bonds: " << num_bonds << std::endl;

    // Allocate memory for bonds
    std::vector<Bond> bonds(num_bonds);  // Vector to hold bond information

    std::vector<double> p2data;

    //Define bond data:
    for (int i = 0; i < num_bonds; i++) {
             imol = i / ( monomers_per_molecule - 2 );
             int remainder = i % ( monomers_per_molecule - 2 );
             iatom = imol * monomers_per_molecule + remainder + 1;
             bonds[i].atom_index = iatom;  // Assign the atom index
     }

    // Read frames
    int frame_count = 0;
    while (read_xtc(xd, natoms, &step, &time, box, x, &prec) == exdrOK) {
	if (frame_count % frame_skip == 0) {
          std::cout << "Time: " << time << " ps" << std::endl;
	  int crystal_count = 0;	// This is the counter to count the bonds that are within a crystalline domain
          for (int i = 0; i < num_bonds; i++) {
		// Assign 3D vector for the bond:
		iatom = bonds[i].atom_index;

		float dx = x[iatom+1][0] - x[iatom-1][0];
		float dy = x[iatom+1][1] - x[iatom-1][1];
		float dz = x[iatom+1][2] - x[iatom-1][2];

     		dx -= round(dx / box[0][0]) * box[0][0];  // Adjust x-coordinate
     		dy -= round(dy / box[1][1]) * box[1][1];  // Adjust y-coordinate
     		dz -= round(dz / box[2][2]) * box[2][2];  // Adjust z-coordinate

		bonds[i].bond_vector[0] = dx;
		bonds[i].bond_vector[1] = dy;
		bonds[i].bond_vector[2] = dz;
		bonds[i].mag = std::sqrt(dx * dx + dy * dy + dz * dz);
		bonds[i].p2 = 0.0;
		bonds[i].neighbors.clear();
	  }
	  for (int i = 0; i < num_bonds; i++) {
		float idx = bonds[i].bond_vector[0];
		float idy = bonds[i].bond_vector[1];
		float idz = bonds[i].bond_vector[2];
		for (int j = 0; j < num_bonds; j++) {
			float jdx = bonds[j].bond_vector[0];
			float jdy = bonds[j].bond_vector[1];
			float jdz = bonds[j].bond_vector[2];
			if ( i < j ) {
				iatom = bonds[i].atom_index;
				jatom = bonds[j].atom_index;
				float dist2 = calculate_distance2_MIC(x[iatom], x[jatom], box);
				if ( dist2 < rcut2 ) {
					bonds[i].neighbors.push_back( j );
					bonds[j].neighbors.push_back( i );
					float costheta = (idx*jdx+idy*jdy+idz*jdz) / (bonds[i].mag * bonds[j].mag);
					bonds[i].p2 += (3.* costheta * costheta - 1.0 ) / 2.0;
					bonds[j].p2 += (3.* costheta * costheta - 1.0 ) / 2.0;
				}
			}
		}
	  }
	  for (int i = 0; i < num_bonds; i++) {
		if (!bonds[i].neighbors.empty()) {
			bonds[i].p2 = bonds[i].p2 / bonds[i].neighbors.size();
		} else {
			bonds[i].p2 = 0;
		}
		p2data.push_back( bonds[i].p2 );
		bonds[i].crystal = false;
		if ( bonds[i].p2 > p2thresh ) {
			crystal_count++;	// If the p2 value for bond i is larger than the threshold value, the bond is labeled as crystal
			bonds[i].crystal = true;
		}
	  }
	std::cout << "Crystalline fraction = " << (float(crystal_count) / num_bonds) << std::endl;
	outputfile << time << " " << (float(crystal_count) / num_bonds) << std::endl;

	// Now we need to perform the clustering analysis:
	// Create a filtered list of only crystalline bonds:
	vector<int> crystal_bond_indices;
	for (int i = 0; i < bonds.size(); ++i) {
	    if (bonds[i].crystal) {
	        crystal_bond_indices.push_back(i);
	    }
	}
	// Compute distances between crystalline bonds only:
	int M = crystal_bond_indices.size();
	vector<vector<int>> adj(M);  // adjacency list for crystalline bonds

	for (int i = 0; i < M; ++i) {
	    iatom = bonds[crystal_bond_indices[i]].atom_index;
	    for (int j = i + 1; j < M; ++j) {
	        jatom = bonds[crystal_bond_indices[j]].atom_index;
	        float dist2 = calculate_distance2_MIC(x[iatom], x[jatom], box);
	        if (dist2 < clustcutoff2) {
	            adj[i].push_back(j);
	            adj[j].push_back(i);
	        }
	    }
	}

	// Cluster Using Iterative DFS:
	
	vector<bool> visited(M, false);
	vector<vector<int>> clusters;

	for (int i = 0; i < M; ++i) {
	    if (!visited[i]) {
	        vector<int> cluster;
	        dfs_iterative(i, adj, visited, cluster);
	        clusters.push_back(cluster);
	    }
	}

	// Pair of (original index, cluster size)
	std::vector<std::pair<int, size_t>> cluster_sizes;
	for (size_t i = 0; i < clusters.size(); ++i) {
	    cluster_sizes.emplace_back(i, clusters[i].size());
	}

	// Sort by size in descending order
	std::sort(cluster_sizes.begin(), cluster_sizes.end(),
	          [](const auto& a, const auto& b) {
	              return a.second > b.second;  // descending
	          });

	std::unordered_map<int, int> old_to_new_id;
	for (size_t new_id = 0; new_id < cluster_sizes.size(); ++new_id) {
	    int old_id = cluster_sizes[new_id].first;
	    old_to_new_id[old_id] = new_id + 1;  // +1 if you want cluster IDs starting from 1
	}

        std::vector<int> atom_cluster(natoms, 0);
	// In this case, each cluster will be a group of indices into the crystal_bond_indices list, so to map back to original center atom indices:
	for (auto& cluster : clusters) {
	    for (int& idx : cluster) {
	        idx = bonds[crystal_bond_indices[idx]].atom_index;  // map back to original center atom index
	    }
	}

	// Printing the Clusters 
	for (size_t i = 0; i < clusters.size(); ++i) {
	    cout << "Cluster " << i + 1 << " (size = " << clusters[i].size() << "): ";
	    for (int idx : clusters[i]) {
	        cout << idx << " ";
		atom_cluster[idx] = old_to_new_id[i];
	    }
	    cout << endl;
	}

	xyz_file << natoms << "\n";
	xyz_file << "Frame " << frame_count << "\n";

	for (int i = 0; i < natoms; ++i) {
	    int cid = atom_cluster[i];
	    std::string label = "c" + std::to_string(cid);  // e.g., c0, c1, ...
	    xyz_file << label << " "
	             << std::fixed << std::setprecision(5)
	             << x[i][0] << " " << x[i][1] << " " << x[i][2] << "\n";
	}

	}
	frame_count++;
    }

    outputfile.close();
    xyz_file.close();

    printHistogram(p2data, binSize);

    // Clean up
    free(x);
    xdrfile_close(xd);

    std::cout << "Processing completed successfully." << std::endl;

    return 0;
}

