import numpy as np
import matplotlib.pyplot as plt

# Set font to serif as a close alternative to Times New Roman
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14  # General font size

# Set larger font sizes specifically for axis labels and tick labels
plt.rcParams["axes.labelsize"] = 18  # Axis label font size
plt.rcParams["xtick.labelsize"] = 16  # X-axis tick label font size
plt.rcParams["ytick.labelsize"] = 16  # Y-axis tick label font size

Npart = 500

# Load simulation data
data = np.loadtxt("thermo_output.txt", skiprows=2)  # Skip header
steps = data[:, 0] / Npart  # MC Steps
TotEng = data[:, 1] / Npart   # Total Energy per particle
Press = data[:, 2]  # System pressure
Acc = data[:, 3]  # Acceptance ratio

# Plot Number of MC steps vs. Energy/N

plt.plot(steps / 10000, TotEng, 'ko', ms=2.5, label=" ")
plt.xlabel("Number of MC steps (" +r'$\times 10^{4}$'+")")
plt.ylabel("Energy per particle (" +r'$\epsilon$'+")")
#plt.legend()
plt.tight_layout()

plt.savefig("energy_vs_time.pdf")  # Save as PDF
plt.close()

# Plot Number of MC steps vs. System pressure

plt.plot(steps / 10000, Press, 'ko', ms=2.5, label=" ")
plt.xlabel("Number of MC steps (" +r'$\times 10^{4}$'+")")
plt.ylabel("Presure (" +r'$\epsilon$'+")")
#plt.legend()
plt.tight_layout()

plt.savefig("pressure_vs_time.pdf")  # Save as PDF
plt.close()

# Plot T vs. Acceptance Ratio

plt.plot(steps / 10000, Acc, 'ko', ms=2.5, label=" ")
plt.xlabel("Number of MC steps ("  +r'$\times 10^{4}$'+")")
plt.ylabel("Acceptance ratio")
#plt.legend()
plt.tight_layout()

plt.savefig("acceptancefrac_vs_time.pdf")  # Save as PDF
plt.close()

#plt.show()

print("Plots saved as energy_vs_time.pdf, pressure_vs_time.pdf, and acceptancefrac_vs_time.pdf")


