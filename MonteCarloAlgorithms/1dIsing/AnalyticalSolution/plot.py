import numpy as np
import matplotlib.pyplot as plt

# Set font to serif as a close alternative to Times New Roman
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14  # General font size

# Set larger font sizes specifically for axis labels and tick labels
plt.rcParams["axes.labelsize"] = 18  # Axis label font size
plt.rcParams["xtick.labelsize"] = 16  # X-axis tick label font size
plt.rcParams["ytick.labelsize"] = 16  # Y-axis tick label font size

# Load simulation data
data = np.loadtxt("results.out", skiprows=1)  # Skip header
T = data[:, 0]    # Temperature (kT/J)
E_N = data[:, 1]  # Energy per spin
Cv_N = data[:, 2] # Heat capacity per spin
Acc = data[:, 3]  # Acceptance fraction

# Define theoretical functions
def theoretical_energy(T):
    return -np.tanh(1.0 / T)  # J = 1, kB = 1

def theoretical_Cv(T):
    return (1.0 / T**2) * (1 / np.cosh(1.0 / T) ** 2)  # J = 1, kB = 1

# Generate theoretical values
T_theory = np.linspace(min(T), max(T), 100)  # Smooth temperature range
E_theory = theoretical_energy(T_theory)
Cv_theory = theoretical_Cv(T_theory)

# Plot T vs. Energy/N
#plt.figure(figsize=(8, 6))

plt.plot(T, E_N, 'ko', ms=12, label="MC simulation data")
plt.plot(T_theory, E_theory, 'k-', label="Exact solution")  # Theoretical curve
plt.xlabel("Temperature ( " + r'$k_{B}T$' + " )")
plt.ylabel("Energy per spin (E/N)")
plt.legend()
plt.tight_layout()

plt.savefig("energy_vs_T.pdf")  # Save as PDF
plt.close()

# Plot T vs. Cv/N
plt.plot(T, Cv_N, 'ko', ms=12, label="MC simulation data")
plt.plot(T_theory, Cv_theory, 'k-', label="Exact solution")  # Theoretical curve
plt.xlabel("Temperature ( " + r'$k_{B}T$' + " )")
plt.ylabel("Heat capacity per spin ( " + r'$C_{V}/N$' + " )")
plt.legend()
plt.tight_layout()

plt.savefig("Cv_vs_T.pdf")  # Save as PDF
plt.close()

# Plot T vs. Acceptance Ratio
plt.plot(T, Acc, 'ko', ms=12, label="Acceptance Fraction")
plt.xlabel("Temperature ( " + r'$k_{B}T$' + " )")
plt.ylabel("Acceptance ratio")
#plt.legend()
plt.tight_layout()

plt.savefig("acceptancefrac_vs_T.pdf")  # Save as PDF
plt.close()

#plt.show()

print("Plots saved as energy_vs_T.pdf, cv_vs_T.pdf, and acceptance_vs_T.pdf")


