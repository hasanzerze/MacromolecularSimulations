import matplotlib.pyplot as plt
import numpy as np

# Set font to serif as a close alternative to Times New Roman
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14  # General font size

# Set larger font sizes specifically for axis labels and tick labels
plt.rcParams["axes.labelsize"] = 18  # Axis label font size
plt.rcParams["xtick.labelsize"] = 16  # X-axis tick label font size
plt.rcParams["ytick.labelsize"] = 16  # Y-axis tick label font size

# Load data from a text file
data = np.loadtxt("ree_vs_n_eps-1.dat", comments="#")

N = data[:, 0]      # Chain lengths
Ree2 = data[:, 1]   # ⟨Ree²⟩
Rg2 = data[:, 2]    # ⟨Rg²⟩

# Expected SARW scaling: R^2 ~ N^{2ν} with ν = 0.75 in 2D
nu = 0.75
N_fit = np.linspace(min(N), max(N), 200)
Ree2_sarw = 0.7764 * N_fit**(2 * nu)    #From Madras & Sokal (1988)
Rg2_sarw = 0.1090* (N_fit**(2 * nu))    #From Madras & Sokal (1988)



plt.figure(figsize=(8, 6))
#plt.ylim(0.5,6.5)
#plt.ylim(0.5,4.8)

# Plot ⟨Ree²⟩
plt.plot(np.log(N), np.log(Ree2), 's', ms=10, label=r'$\langle R_{ee}^2 \rangle$', color='black')

# Plot ⟨Rg²⟩
plt.plot(np.log(N), np.log(Rg2), '^', ms=10, label=r'$\langle R_g^2 \rangle$', color='black')

# Plot expected SARW power law as a reference
#plt.plot(np.log(N_fit), np.log(Ree2_sarw), 'k-', label=r'$\langle R_{ee}^2 \rangle$'+", Madras & Sokal (1988)")
#plt.plot(np.log(N_fit), np.log(Rg2_sarw), 'k--', label=r'$\langle R_g^2 \rangle$'+", Madras & Sokal (1988)")
plt.plot(np.log(N_fit), np.log(N_fit), 'k-', label=r'$\langle R^2 \rangle$' + " ~ N")

plt.xlabel("$\log(N)$")
plt.ylabel(r'$\log\langle R^2 \rangle$')
plt.legend()
plt.tight_layout()

plt.savefig("scaling_plot_eps-1.png", dpi=300)
#plt.show()

